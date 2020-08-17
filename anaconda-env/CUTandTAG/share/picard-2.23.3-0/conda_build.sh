if [ -z ${CONDA_BUILD+x} ]; then
    source /opt/conda/conda-bld/picard_1595279665491/work/build_env_setup.sh
fi
#!/bin/bash

TGT="$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM"
[ -d "$TGT" ] || mkdir -p "$TGT"
[ -d "${PREFIX}/bin" ] || mkdir -p "${PREFIX}/bin"

cd "${SRC_DIR}"
# Do not install Linux specific x86-acceleration libraries
if [ "$(uname)" == "Darwin" ]; then
    rm -f libIntel*.so
fi
cp -rvp . "${TGT}"

cp $RECIPE_DIR/picard.sh $TGT/picard
ln -s $TGT/picard $PREFIX/bin
chmod 0755 "${PREFIX}/bin/picard"
