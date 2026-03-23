from sgffp import SgffReader, SgffObject, SgffSegment, SgffFeature
from sgffp.models.history import SgffHistoryNode, SgffHistoryTreeNode
from Bio.Seq import Seq
import re
from pydna.assembly2 import gibson_assembly, pcr_assembly, restriction_ligation_assembly, gateway_assembly, in_fusion_assembly
from pydna.primer import Primer
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from pydna.dseqrecord import Dseqrecord
import glob
import os
from pydna.opencloning_models import CloningStrategy, AssemblyFragment, Source
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Restriction import RestrictionBatch
from pydna.parsers import parse_snapgene
from pydna.utils import flatten


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

    seq = Seq(node.sequence)
    circular = tree_node.circular if tree_node else False
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
        record=str(seq),
        id=name,
        name=name,
        description=name,
        features=features,
        annotations=annotations,
        circular=circular,
    )
    return record



def filter_assembly_fragments_that_are_sequences(input_value: list[AssemblyFragment]) -> list[AssemblyFragment]:
    return [fragment for fragment in input_value if isinstance(fragment.sequence, Dseqrecord)]

def source_from_tree_node(
    product: Dseqrecord, node: SgffHistoryTreeNode, sgff_object: SgffObject
) -> tuple[Source|None|int, list[SgffHistoryNode]]:
    input_sequences = [
        history_node_to_dseqrecord(sgff_object, child.id) for child in node.children
    ]
    if len(input_sequences) == 0:
        return None, []
    expected_seguid = product.seq.seguid()
    if node.operation == "gibsonAssembly":
        products = gibson_assembly(input_sequences, limit=15)
    elif node.operation == "amplifyFragment":
        primers = [Primer(oligo.sequence, name=oligo.name) for oligo in node.oligos]
        products = pcr_assembly(input_sequences[0], *primers, limit=12)
    elif node.operation == "changeStrandedness":
        return -1, None
    elif node.operation == "insertFragment" or node.operation == "goldenGateAssembly":
        enzyme_names = set(flatten([input_summary.enzyme_names for input_summary in node.input_summaries]))
        if all(enz_name in rest_dict.keys() for enz_name in enzyme_names):
            rb = RestrictionBatch(first=[e for e in enzyme_names])
            products = restriction_ligation_assembly(input_sequences, rb)
        else:
            raise ValueError(f"Unknown enzymes: {enzyme_names}")
    elif node.operation == "gatewayLRCloning":
        products = gateway_assembly(input_sequences, "LR")
    elif node.operation == "gatewayBPCloning":
        products = gateway_assembly(input_sequences, "BP")
    elif node.operation == "inFusionCloning":
        products = in_fusion_assembly(input_sequences, limit=15)
    else:
        raise ValueError(f"Unknown operation: {node.operation}")

    correct_product = next(
        (p for p in products if p.seq.seguid() == expected_seguid), None
    )
    if correct_product is None:
        raise ValueError(f"No product found for expected SEGUID {expected_seguid}")

    # Return the children nodes in the same order as the inputs
    out_nodes = [None] * len(input_sequences)
    for input_value in filter_assembly_fragments_that_are_sequences(correct_product.source.input):
        idx = next((i for i, seq in enumerate(input_sequences) if seq == input_value.sequence), None)
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
    for input_value in filter_assembly_fragments_that_are_sequences(source.input):
        node = out_nodes.pop(0)
        parse_history(input_value.sequence, node, sgff_object)


# --- Test it ---
for file in glob.glob("data/*.dna"):
    print(file)
    root_record = parse_snapgene(file)[0]
    root_record.name = os.path.basename(file)
    sgff_object = SgffReader.from_file(file)
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
