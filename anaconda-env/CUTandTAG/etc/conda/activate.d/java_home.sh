export CONDA_BACKUP_JAVA_HOME="${JAVA_HOME}"
export JAVA_HOME="${CONDA_PREFIX}"

export CONDA_BACKUP_JAVA_LD_LIBRARY_PATH="${JAVA_LD_LIBRARY_PATH}"

if [ $(uname) = Darwin ]; then
  export JAVA_LD_LIBRARY_PATH="${JAVA_HOME}"/jre/lib/server
else
  if [ $(uname -m) = x86_64 ]; then
    export JAVA_LD_LIBRARY_PATH="${JAVA_HOME}"/jre/lib/amd64/server
  else
    export JAVA_LD_LIBRARY_PATH="${JAVA_HOME}"/jre/lib/i386/server
  fi
fi
