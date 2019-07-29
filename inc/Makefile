default: amgcl taskflow

.PHONY: amgcl taskflow

amgcl:
	rm -fr amgcl
	-svn export https://github.com/ddemidov/amgcl.git/trunk/amgcl
	@if [ ! -d "amgcl" ]; then  unzip amgcl.zip; echo "amgcl svn export failed: unzipped local copy."; fi

taskflow:
	rm -fr taskflow
	-svn export https://github.com/cpp-taskflow/cpp-taskflow.git/trunk/taskflow
	@if [ ! -d "taskflow" ]; then  unzip taskflow.zip; echo "taskflow svn export failed: unzipped local copy."; fi

