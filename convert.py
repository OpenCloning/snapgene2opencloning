from sgffp import SgffReader, SgffObject, SgffSegment, SgffFeature
from sgffp.models.history import (
    SgffHistoryNode,
    SgffHistoryTreeNode,
    SgffInputSummary,
    SgffHistoryOligo,
)
from pydna.dseq import Dseq
import re
from pydna.assembly2 import (
    gibson_assembly,
    pcr_assembly,
    restriction_ligation_assembly,
    gateway_assembly,
    in_fusion_assembly,
    ligation_assembly,
    fusion_pcr_assembly,
)
from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization
from pydna.primer import Primer
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from pydna.dseqrecord import Dseqrecord
import glob
import os
from pydna.opencloning_models import (
    CloningStrategy,
    AssemblyFragment,
    Source,
)
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Restriction import RestrictionBatch
from pydna.parsers import parse_snapgene
from pydna.utils import flatten
import itertools
import warnings


STRAND_MAP = {"+": 1, "-": -1, ".": 0, "=": 0}


def _segments_to_location(
    segments: list[SgffSegment], strand_int: int, seq_len: int, circular: bool
) -> SimpleLocation | CompoundLocation | None:
    """Convert SgffSegment list to a Biopython location."""
    locations = []
    for seg in segments:
        if seg.type == "gap":
            continue
        # Handle wrap-around for circular sequences
        if circular and seg.start > seg.end:
            locations.append(SimpleLocation(seg.start, seq_len, strand=strand_int))
            locations.append(SimpleLocation(0, seg.end, strand=strand_int))
        else:
            locations.append(SimpleLocation(seg.start, seg.end, strand=strand_int))

    if not locations:
        return None
    if len(locations) == 1:
        return locations[0]
    # For reverse strand, reverse order
    if strand_int == -1:
        locations = locations[::-1]
    return CompoundLocation(locations)


def _feature_to_seqfeature(
    feature: SgffFeature, seq_len: int, circular: bool
) -> SeqFeature | None:
    """Convert an SgffFeature to a Biopython SeqFeature."""
    strand_int = STRAND_MAP.get(feature.strand, 0)
    location = _segments_to_location(feature.segments, strand_int, seq_len, circular)
    if location is None:
        return None

    # Convert qualifiers: Biopython expects lists as values
    qualifiers = {
        k: [v] if not isinstance(v, list) else v for k, v in feature.qualifiers.items()
    }
    if feature.name and "label" not in qualifiers:
        qualifiers["label"] = [feature.name]

    return SeqFeature(location=location, type=feature.type, qualifiers=qualifiers)


def history_node_to_dseqrecord(sgff_object: SgffObject, node_id: str) -> Dseqrecord:
    """Convert a history node to a Dseqrecord.

    Sequence comes from the history node, metadata from the tree node.
    """
    node: SgffHistoryNode = sgff_object.history.nodes[node_id]
    tree_node = sgff_object.history.get_tree_node(node_id)

    circular = tree_node.circular if tree_node else False
    seq_props = node.properties.get("AdditionalSequenceProperties")
    if circular:
        seq = Dseq(node.sequence, circular=True)
    elif (
        seq_props is not None
        and "UpstreamStickiness" in seq_props
        and "DownstreamStickiness" in seq_props
    ):
        left_ovhg = -int(seq_props.get("UpstreamStickiness"))
        right_ovhg = -int(seq_props.get("DownstreamStickiness"))
        seq = Dseq.from_full_sequence_and_overhangs(
            node.sequence, left_ovhg, right_ovhg
        )
    else:
        seq = Dseq(node.sequence)
    seq_len = node.length
    name = tree_node.name if tree_node else f"node_{node_id}"
    name = re.sub(r"\s+", "_", name)  # Replace whitespace with underscores

    # Convert features from the node's content
    features = []
    for feat in node.features:
        sf = _feature_to_seqfeature(feat, seq_len, circular)
        if sf is not None:
            features.append(sf)

    annotations = {}
    annotations["topology"] = "circular" if circular else "linear"
    if tree_node:
        annotations["molecule_type"] = "DNA"
        if tree_node.strandedness:
            annotations["molecule_type"] = f"{tree_node.strandedness}-stranded DNA"

    record = Dseqrecord(
        record=seq,
        id=name,
        name=name,
        description=name,
        features=features,
        annotations=annotations,
    )
    return record


def filter_assembly_fragments_that_are_sequences(input_value: list[AssemblyFragment]) -> list[AssemblyFragment]:
    return [fragment for fragment in input_value if isinstance(fragment.sequence, Dseqrecord)]


def get_enzyme_batch_from_input_summaries(
    input_summaries: list[SgffInputSummary],
) -> RestrictionBatch:
    enzyme_names = set(
        flatten([input_summary.enzyme_names for input_summary in input_summaries])
    )
    # Sometimes enzymes come with < > around, we remove them
    enzyme_names = set(
        enz_name.replace("<", "").replace(">", "") for enz_name in enzyme_names
    )
    # Sometimes they have Start or End tags, we remove them
    enzyme_names = enzyme_names.difference({"Start", "End"})
    if all(enz_name in rest_dict.keys() for enz_name in enzyme_names):
        return RestrictionBatch(first=[e for e in enzyme_names])
    else:
        raise ValueError(f"Unknown enzymes: {enzyme_names}")


def get_sequence_inputs(source: Source) -> list[Dseqrecord]:
    """Auxiliary function to get the most ancestral sequences used as inputs. These will tipically be the immediate inputs,
    but in a case where a restriction-ligation is turned into restriction, then ligation, we have to go up the tree to find the most ancestral sequences.
    """
    out_value = list()
    for input_value in filter_assembly_fragments_that_are_sequences(source.input):
        if (
            input_value.sequence.source is None
            or len(input_value.sequence.source.input) == 0
        ):
            out_value.append(input_value.sequence)
        else:
            out_value.extend(get_sequence_inputs(input_value.sequence.source))
    return out_value


def parseOligos(oligos: list[SgffHistoryOligo]) -> list[Primer]:
    return [
        Primer(oligo.sequence, name=oligo.name or f"oligo_{i+1}")
        for i, oligo in enumerate(oligos)
    ]


def source_from_tree_node(
    expected_product: Dseqrecord, node: SgffHistoryTreeNode, sgff_object: SgffObject
) -> tuple[Source | None | int, list[SgffHistoryNode]]:
    input_sequences = [
        history_node_to_dseqrecord(sgff_object, child.id) for child in node.children
    ]

    expected_seguid = expected_product.seq.seguid()
    expected_dseq_and_rc = (
        expected_product.seq,
        expected_product.seq.reverse_complement(),
    )

    def find_expected_product(products: list[Dseqrecord]) -> Dseqrecord | None:
        if expected_product.circular:
            return next(
                (p for p in products if p.seq.seguid() == expected_seguid), None
            )
        else:
            return next((p for p in products if p.seq in expected_dseq_and_rc), None)

    print(">>", node.operation)

    if node.operation == "gibsonAssembly":
        products = gibson_assembly(input_sequences, limit=15)
    elif node.operation == "amplifyFragment":
        primers = parseOligos(node.oligos)
        products = pcr_assembly(input_sequences[0], *primers, limit=12)
    elif node.operation == "primerDirectedMutagenesis":
        fwd_primer, *_ = parseOligos(node.oligos)
        rvs_primer = Primer(
            fwd_primer.seq.reverse_complement(), name=f"rvs_{fwd_primer.name}"
        )
        pcr_products = pcr_assembly(
            input_sequences[0], fwd_primer, rvs_primer, limit=10
        )
        # They bundle also the fusion pcr on the same step
        products = list()
        for pcr_product in pcr_products:
            pcr_product.name = f"mutagenesis_pcr_product"
            products.extend(fusion_pcr_assembly([pcr_product], limit=6))
    elif node.operation in ["changeStrandedness", "editDNAEnds", "changeMethylation"]:
        return -1, None
    elif node.operation == "changeTopology":
        if expected_product.circular:
            input_sequences = [expected_product[: len(expected_product)]]
            products = ligation_assembly(input_sequences, True)
            if len(products) == 0:
                warnings.warn(f"Stopped at change topology operation")
                return None, []
        else:
            warnings.warn(f"Stopped at change topology operation")
            return None, []
    elif node.operation in [
        "insertFragment",
        "goldenGateAssembly",
        "insertFragments",
        "ligateFragments",
    ]:

        rb = get_enzyme_batch_from_input_summaries(node.input_summaries)
        products = restriction_ligation_assembly(input_sequences, rb)
        if find_expected_product(products) is None:
            # Try a simple ligation if the restriction ligation failed
            products = ligation_assembly(input_sequences)
        if find_expected_product(products) is None:
            # Try restriction, then ligation if the simple ligation failed
            digestion_products = list()
            for input_sequence, input_summary in zip(
                input_sequences, node.input_summaries
            ):
                rb = get_enzyme_batch_from_input_summaries([input_summary])
                digestion_products.append(input_sequence.cut(rb))

            possible_combinations = itertools.product(*digestion_products)
            for combination in possible_combinations:
                products = ligation_assembly(combination)
                if find_expected_product(products) is not None:
                    break

    elif node.operation == "gatewayLRCloning":
        products = gateway_assembly(input_sequences, "LR")
    elif node.operation == "gatewayBPCloning":
        products = gateway_assembly(input_sequences, "BP")
    elif node.operation == "inFusionCloning":
        products = in_fusion_assembly(input_sequences, limit=10)
    elif node.operation == "hifiAssembly":
        products = gibson_assembly(input_sequences, limit=10)
    elif node.operation == "overlapFragments":
        products = fusion_pcr_assembly(input_sequences, limit=10)
    elif node.operation == "annealOligos":
        primers = parseOligos(node.oligos)
        products = oligonucleotide_hybridization(*primers, 10)
    elif node.operation == "invalid":
        return None, []
    else:
        raise ValueError(f"Unknown operation: {node.operation}")

    correct_product = find_expected_product(products)

    if correct_product is None:
        raise ValueError(f"No product found for expected SEGUID {expected_seguid}")

    # Return the children nodes in the same order as the inputs
    out_nodes = [None] * len(input_sequences)
    for input_value in get_sequence_inputs(correct_product.source):
        idx = next(
            (i for i, seq in enumerate(input_sequences) if seq is input_value),
            None,
        )
        out_nodes[idx] = node.children[idx]

    return correct_product.source, out_nodes


def parse_history(root_record: Dseqrecord, root_node: SgffHistoryTreeNode, sgff_object: SgffObject) -> Dseqrecord:
    """Parse the history of a Dseqrecord, and edit it in place."""
    repeat = True
    while repeat:
        source, out_nodes = source_from_tree_node(root_record, root_node, sgff_object)
        repeat = source == -1
        if repeat:
            root_node = root_node.children[0]

    root_record.source = source
    if source is None:
        return
    for input_value in get_sequence_inputs(source):
        node = out_nodes.pop(0)
        parse_history(input_value, node, sgff_object)


# --- Test it ---
for file in glob.glob("data/*.dna"):
    # for file in glob.glob("data/linear_ligation_overhangs2.dna"):
    print(file)

    root_record = parse_snapgene(file)[0]
    root_record.name = os.path.basename(file)
    sgff_object = SgffReader.from_file(file)
    seq_props = sgff_object.properties.get("AdditionalSequenceProperties")
    if (
        not root_record.circular
        and seq_props is not None
        and "UpstreamStickiness" in seq_props
        and "DownstreamStickiness" in seq_props
    ):
        left_ovhg = -int(seq_props.get("UpstreamStickiness"))
        right_ovhg = -int(seq_props.get("DownstreamStickiness"))
        root_record.seq = Dseq.from_full_sequence_and_overhangs(
            str(root_record.seq), left_ovhg, right_ovhg
        )

    if not sgff_object.has_history:
        print("No history found, skipping")
        continue
    parse_history(root_record, sgff_object.history.tree.root, sgff_object)
    root_record = root_record.normalize_history()
    root_record.validate_history()
    cs = CloningStrategy.from_dseqrecords([root_record])
    file_name = os.path.basename(file)
    with open(f"output/{file_name.replace('.dna', '.json')}", "w") as f:
        f.write(cs.model_dump_json(indent=2))
