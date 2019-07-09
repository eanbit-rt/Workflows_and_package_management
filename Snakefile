rule complement:
  input:
    "data/{file}.txt"
  output:
    "complement_data/{file}.txt"
  shell:
    "cat {input}|tr atcg tagc > {output}"

rule reverse:
  input:
    "complement_data/{file}.txt"
  output:
    "reverse_data/{file}.md"
  shell:
    "cat {input}|rev >{output}"
