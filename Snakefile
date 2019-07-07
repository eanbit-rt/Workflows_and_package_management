rule complement:
  input:
    "data/text.txt"
  output:
    "complement_data/text.txt"
  shell:
    "cat {input}|tr atcg tagc > {output}"
