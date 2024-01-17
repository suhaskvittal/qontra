#!/bin/sh

CODE=$1
ROUNDS=$2

cd Release

TANNER_GRAPH_FILE=../data/tanner/${CODE}.txt
ARCH_FOLDER_PREFIX=../data/protean/${CODE}

# Version 1: no modifications.
./protean --tanner $TANNER_GRAPH_FILE \
            --out ${ARCH_FOLDER_PREFIX}/v1 \
            --passes "ral.rlb.rcr" \
            --s-rounds $ROUNDS \
            --render ${ARCH_FOLDER_PREFIX}/v1/render
# Version 2.1: add flags
./protean --tanner $TANNER_GRAPH_FILE \
            --out ${ARCH_FOLDER_PREFIX}/v2.1 \
            --passes "fla.ral.rlb.rcr" \
            --s-rounds $ROUNDS \
            --render ${ARCH_FOLDER_PREFIX}/v2.1/render
# Version 2.2: join identical supports + add flags
./protean --tanner $TANNER_GRAPH_FILE \
            --out ${ARCH_FOLDER_PREFIX}/v2.2 \
            --passes "jid.ral.fla.ral.rlb.rcr" \
            --s-rounds $ROUNDS \
            --render ${ARCH_FOLDER_PREFIX}/v2.2/render
# Version 3.1: v2.1 + other optimizations.
./protean --tanner $TANNER_GRAPH_FILE \
            --out ${ARCH_FOLDER_PREFIX}/v3.1 \
            --passes "fla.ral(con.ral.prx.ral)+rlb.rcr" \
            --s-rounds $ROUNDS \
            --render ${ARCH_FOLDER_PREFIX}/v3.1/render
# Version 3.2: v2.2 + other optimizations.
./protean --tanner $TANNER_GRAPH_FILE \
            --out ${ARCH_FOLDER_PREFIX}/v3.2 \
            --passes "jid.ral.fla.ral(con.ral.prx.ral)+rlb.rcr" \
            --s-rounds $ROUNDS \
            --render ${ARCH_FOLDER_PREFIX}/v3.2/render
