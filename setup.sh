for i in memod-s  Snakefile  scripts  envs
do
        cp -r $i $CONDA_PREFIX/bin
        chmod -R +x $CONDA_PREFIX/bin/$i
done
