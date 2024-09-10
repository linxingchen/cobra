from cobra import *
from cobra import _get_subset_trunks

query_set = frozenset({"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"})


def test_get_sub_trunks_one():
    # fmt: off
    contig2assembly = {
        "A": {"A", "B"},
        "B": {"A", "B"},
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
    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset(),
        contig_link_no_pe=frozenset(),
    )
    print(assembly_reason)
    assert len(assembly_reason) == 1
    assert list(assembly_reason.values())[0] == AssemblyReason(
        groupid=0,
        judgement="standalone",
        represent_seqs=[],
        dup_queries=frozenset({"A", "B"}),
    )


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

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset(),
        contig_link_no_pe=frozenset(),
    )
    print(assembly_reason)
    assert assembly_reason == {
        "B": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=assembly_reason["A"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "A": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=assembly_reason["B"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "C": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=assembly_reason["C"].represent_seqs,
            dup_queries=frozenset(),
        ),
    }
    assert len(assembly_reason["A"].represent_seqs) == 1
    assert set(assembly_reason["C"].represent_seqs) == {"A", "B"}


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

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset("D"),
        contig_link_no_pe=frozenset(),
    )
    print(assembly_reason)
    assert assembly_reason == {
        "D": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=["D"],
            dup_queries=frozenset(),
        )
    }

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset(),
        contig_link_no_pe=frozenset("C"),
    )
    print(assembly_reason)
    # FIXME: this can be a bug: we expect at least A--B can be kept
    assert assembly_reason == {
        "D": AssemblyReason(
            groupid=0,
            judgement="nolink_query",
            represent_seqs=["D"],
            dup_queries=frozenset(),
        )
    }

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset("B"),
        contig_link_no_pe=frozenset(),
    )
    print(assembly_reason)
    assert assembly_reason == {
        "D": AssemblyReason(
            groupid=0,
            judgement="circular_in_sub",
            represent_seqs=["D"],
            dup_queries=frozenset(),
        ),
        "C": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=["C"],
            dup_queries=frozenset(),
        ),
        "B": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=["A", "B"],
            dup_queries=frozenset(),
        ),
    }


def test_get_sub_trunks_disjoint_i():
    # fmt: off
    contig2assembly = {
        "A": {"A",              },
        "B": {"A", "B"          },
        "C": {          "C", "D"},
        "D": {     "B", "C", "D"},
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
        "C": SubsetChunk(standalong_subs={}, frags=["C"]),
        "A": SubsetChunk(standalong_subs={}, frags=["A"]),
    }
    assert _unextendable == {"D", "B"}
    assert _failed_reason == {
        "D": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=_failed_reason["B"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "B": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=_failed_reason["D"].represent_seqs,
            dup_queries=frozenset(),
        ),
    }


def test_get_sub_trunks_disjoint_sub():
    # fmt: off
    contig2assembly = {
        "A": {"A",                   },
        "B": {"A", "B"               },
        "C": {          "C", "D"     },
        "D": {"A", "B", "C", "D"     },
        "E": {     "B", "C", "D", "E"},
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
        "A": SubsetChunk(standalong_subs={}, frags=["A"]),
        "C": SubsetChunk(standalong_subs={}, frags=["C"]),
    }
    assert len(_unextendable) == 3
    assert _unextendable < {"B", "C", "D", "E"}
    assert _failed_reason == {
        "E": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=_failed_reason["E"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "D": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=_failed_reason["D"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "B": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=["B"],
            dup_queries=frozenset(),
        ),
    }


def test_get_sub_trunks_8():
    # fmt: off
    contig2assembly = {
        "A": {"A",         },
        "B": {"A", "B"     },
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
    print(groups2ext_query)

    subset_trunks, _unextendable, _failed_reason = _get_subset_trunks(
        0, groups2ext_query[0]
    )
    print(f"{subset_trunks=}", f"{_unextendable=}", f"{_failed_reason=}", sep="\n")
    assert subset_trunks == {
        "C": SubsetChunk(standalong_subs={}, frags=["A", "B", "C"])
    }
    assert not set(_unextendable)
    assert not set(_failed_reason)

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset({"B", "C"}),
        contig_link_no_pe=frozenset(),
    )
    print(assembly_reason)
    assert assembly_reason == {
        "C": AssemblyReason(
            groupid=0,
            judgement="circular_8_tight",
            represent_seqs=["C"],
            dup_queries=frozenset(),
        ),
        "B": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=["A", "B"],
            dup_queries=frozenset(),
        ),
    }


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

    assembly_reason = get_assembly2reason(
        0,
        groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset("D"),
        contig_link_no_pe=frozenset(),
    )
    print(assembly_reason)
    assert assembly_reason == {
        "B": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=assembly_reason["A"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "A": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=assembly_reason["B"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "D": AssemblyReason(
            groupid=0,
            judgement="conflict_query",
            represent_seqs=assembly_reason["D"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "C": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=["C"],
            dup_queries=frozenset(),
        ),
    }


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
