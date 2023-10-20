import libVEQ

libVEQ.run(
    [
        "dummy",
        "-dg",
        "graph/data/COLLAB.gfu",
        "-qg",
        "graph/query/COLLAB/randomwalk/8/q30.gfu",
    ],
)
