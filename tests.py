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
    # FIXME: 这里其实应该 OK? 因为 A B 都是 C 的子集
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
        "C": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=assembly_reason["C"].represent_seqs,
            dup_queries=frozenset(),
        ),
    }
    assert set(assembly_reason["C"].represent_seqs) == {"A", "B"}


def test_get_sub_trunks_subs():
    # fmt: off
    contig2assembly = {
        "A": {"A"               },
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
            represent_seqs=assembly_reason["D"].represent_seqs,
            dup_queries=frozenset(),
        )
    }
    assert set(assembly_reason["D"].represent_seqs) == {"A", "B", "C"}

    # 实际上不可能存在
    assembly_reason = get_assembly2reason(
        groupi=0,
        group=groups2ext_query[0],
        contig2assembly=contig2assembly,
        path_circular_potential=frozenset(),
        contig_link_no_pe=frozenset("C"),
    )
    print(assembly_reason)
    # FIXME: this can be a bug: we expect at least A--B can be kept
    assert assembly_reason["D"].judgement == "nolink_query"
    assert assembly_reason["B"].judgement == "longest"
    assert assembly_reason["B"] == AssemblyReason(
        groupid=0, judgement="longest", represent_seqs=["A"], dup_queries=frozenset()
    )

    # FIXME: 应该都失败? 或者丢掉圈, 仅保留臂
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
            judgement="circular_6_conflict",
            represent_seqs=assembly_reason["D"].represent_seqs,
            dup_queries=frozenset(),
        ),
        "C": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=[],
            dup_queries=frozenset(),
        ),
        "A": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=[],
            dup_queries=frozenset(),
        ),
    }
    assert set(assembly_reason["D"].represent_seqs) == {"B", "D"}


def test_get_sub_trunks_disjoint_i():
    # fmt: off
    contig2assembly = {
        "A": {"A"               },
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
        "A": {"A"                    },
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
        "A": {"A"          },
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
            represent_seqs=["A"],
            dup_queries=frozenset(),
        ),
    }


def test_get_sub_trunks_sub_only():
    # FIXME: 这里应该保留 D, ABC 都是子集
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
        "D": AssemblyReason(
            groupid=0,
            judgement="longest",
            represent_seqs=assembly_reason["D"].represent_seqs,
            dup_queries=frozenset(),
        ),
    }
    assert set(assembly_reason["D"].represent_seqs) == {"A", "B", "C"}


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
    (rep_conflit_seq,) = {
        "AcMG_9925",
        "AcMG_5518",
        "AcMG_12242",
    } & assembly_reason.keys()
    assert assembly_reason == {
        rep_conflit_seq: AssemblyReason(
            groupid=0,
            judgement="circular_6_conflict",
            represent_seqs=assembly_reason[rep_conflit_seq].represent_seqs,
            dup_queries=frozenset({"AcMG_9925", "AcMG_5518", "AcMG_12242"}),
        ),
        "AcMG_3793": AssemblyReason(
            groupid=0, judgement="longest", represent_seqs=[], dup_queries=frozenset()
        ),
    }
    assert set(assembly_reason[rep_conflit_seq].represent_seqs) == {
        rep_conflit_seq,
        "AcMG_755",
    }


def _test_local(
    groups2ext_query: dict[int, dict[str, GroupAssemblyIndex]],
    contig2join: dict[str, list[str]],
    contig2cov: dict[str, float],
    contig_pe_links: frozenset[tuple[str, str]],
):
    for groupi, group in groups2ext_query.items():
        if len(group) == 1:
            continue
        contigs: set[str] = set(group)
        contig_links: set[tuple[str, str]] = set()
        for contig in group:
            for end in (f"{contig}_L", f"{contig}_R"):
                last_end = end
                for next_end in contig2join.get(end, []):
                    if next_end.endswith("rc"):
                        next_end = next_end[:-2]
                    contig_links.add((last_end, next_end))
                    last_end = end2end2(next_end)
                    contigs.add(end2contig(next_end))
        gv_names: list[str] = []
        gv_covs: list[str] = []
        gv_links: list[str] = []
        for contig in contigs:
            gv_names.append(f'{contig}_L -> {contig}_R [label="{contig}"; color=blue]')
            gv_covs.append(
                f'{contig}_R -> {contig}_L [label="{contig2cov[contig]}"; color=blue]'
            )
        for pair in contig_links:
            color = ""
            if pair not in contig_pe_links:
                color = "[color=gray]"
            gv_links.append(f"{pair[0]} -> {pair[1]} {color}")
        gv_str = (
            f"digraph group_{groupi} "
            + "{node[shape=box, style=rounded];"
            + (";".join(gv_names) + ";")
            + (";".join(gv_covs) + ";")
            + (";".join(gv_links) + ";")
            + "}"
        )
        print(gv_str)
