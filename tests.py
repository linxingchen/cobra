from cobra import *
from cobra import _get_subset_trunks

query_set = frozenset({"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"})


def test_get_sub_trunks_conflict():
    # fmt: off
    contig2assembly = {
        "A": {"A", "B"     },
        "B": {     "B", "C"},
        "C": {"A", "B", "C"},
    }
    # fmt: on
    groups2ext_query = dict(
        enumerate(
            sorted(
                query2groups(contig2assembly, query_set=query_set),
                key=lambda d: sorted(d),
            )
        )
    )
    subset_trunks, _unextendable, _failed_reason = _get_subset_trunks(
        0, groups2ext_query[0]
    )
    print(f"{subset_trunks=}", f"{_unextendable=}", f"{_failed_reason=}", sep="\n")
    print()
    assert subset_trunks == {}
    assert set(_unextendable) == {"A", "B", "C"}
    assert set(_failed_reason) == {"A", "B", "C"}


def test_get_sub_trunks_subs():
    # fmt: off
    contig2assembly = {
        "A": {"A",              },
        "B": {"A", "B"          },
        "C": {          "C", "D"},
        "D": {"A", "B", "C", "D"},
    }
    # fmt: on
    groups2ext_query = dict(
        enumerate(
            sorted(
                query2groups(contig2assembly, query_set=query_set),
                key=lambda d: sorted(d),
            )
        )
    )
    print(groups2ext_query)

    subset_trunks, _unextendable, _failed_reason = _get_subset_trunks(
        0, groups2ext_query[0]
    )
    print(f"{subset_trunks=}", f"{_unextendable=}", f"{_failed_reason=}", sep="\n")
    assert subset_trunks == {
        "D": SubsetChunk(
            standalong_subs={
                "B": SubsetChunk(standalong_subs={}, frags=["A", "B"]),
                "C": SubsetChunk(standalong_subs={}, frags=["C"]),
            },
            frags=["D"],
        )
    }
    assert not set(_unextendable)
    assert not set(_failed_reason)


def test_get_sub_trunks_sub_only():
    # fmt: off
    contig2assembly = {
        "A": {"A", "B"          },
        "B": {     "B", "C", "D"},
        "C": {          "C", "D"},
        "D": {"A", "B", "C", "D"},
    }
    # fmt: on
    groups2ext_query = dict(
        enumerate(
            sorted(
                query2groups(contig2assembly, query_set=query_set),
                key=lambda d: sorted(d),
            )
        )
    )
    print(groups2ext_query)

    subset_trunks, _unextendable, _failed_reason = _get_subset_trunks(
        0, groups2ext_query[0]
    )
    print(f"{subset_trunks=}", f"{_unextendable=}", f"{_failed_reason=}", sep="\n")
    assert subset_trunks == {"C": SubsetChunk(standalong_subs={}, frags=["C"])}
    assert set(_unextendable) == {"A", "B", "D"}
    assert set(_failed_reason) == {"A", "B", "D"}


def test_get_sub_trunks_real():
    # fmt: off
    query_set = frozenset((           "AcMG_5518", "AcMG_12242", "AcMG_3793", "AcMG_755", "AcMG_9925", "AcMG_9528"))
    contig2assembly = {
        "AcMG_3793" : {"AcMG_925032", "AcMG_5518", "AcMG_12242", "AcMG_3793"},
        "AcMG_755"  : {"AcMG_925032", "AcMG_5518", "AcMG_12242", "AcMG_3793", "AcMG_755", "AcMG_9925"},
        "AcMG_12242": {"AcMG_925032"} | query_set,
        "AcMG_5518" : {"AcMG_925032"} | query_set,
        "AcMG_9925" : {"AcMG_925032"} | query_set,
    }
    # fmt: on
    groups2ext_query = dict(
        enumerate(
            sorted(
                query2groups(contig2assembly, query_set=query_set),
                key=lambda d: sorted(d),
            )
        )
    )
    assert set(groups2ext_query[0].values()) == set(
        {
            "AcMG_3793": GroupAssemblyIndex(
                special=frozenset({"AcMG_12242"}), dup_queries=frozenset()
            ),
            "AcMG_755": GroupAssemblyIndex(
                special=frozenset({"AcMG_755", "AcMG_12242"}), dup_queries=frozenset()
            ),
            "AcMG_5518": GroupAssemblyIndex(
                special=frozenset({"AcMG_755", "AcMG_9528", "AcMG_12242"}),
                dup_queries=frozenset({"AcMG_12242", "AcMG_5518", "AcMG_9925"}),
            ),
        }.values()
    )
    print(groups2ext_query)

    subset_trunks, _unextendable, _failed_reason = _get_subset_trunks(
        0, groups2ext_query[0]
    )
    print(f"{subset_trunks=}", f"{_unextendable=}", f"{_failed_reason=}", sep="\n")
    assert len(subset_trunks) == 1
    (subset,) = subset_trunks.values()
    assert subset.standalong_subs == {}
    assert {"AcMG_3793", "AcMG_755"} < set(subset.frags) and len(subset.frags) == 3
    assert _unextendable == set()
    assert _failed_reason == {}

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset({"AcMG_755"}),
        contig_link_no_pe=frozenset({"AcMG_9528"}),
    )
    print(assembly_reason)
    assert assembly_reason == {
        "AcMG_755": AssemblyReason(
            groupid=0,
            judgement="circular_6_conflict",
            represent_seqs=assembly_reason["AcMG_755"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "AcMG_3793": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=["AcMG_3793"],
            dup_queries=frozenset(),
        ),
    }
    assert len(assembly_reason["AcMG_755"].represent_seqs) == 2
    assert assembly_reason["AcMG_755"].represent_seqs[0] == "AcMG_755"


test_get_sub_trunks_conflict()
test_get_sub_trunks_subs()
test_get_sub_trunks_sub_only()
test_get_sub_trunks_real()
