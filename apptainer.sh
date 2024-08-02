apptainer run \
  --bind $PWD:/data \
  docker://openjournals/inara \
  -o pdf,crossref \
  /data/paper.md
