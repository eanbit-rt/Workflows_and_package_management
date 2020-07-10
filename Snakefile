rule complement:
  input:
    "data/{file}"
  output:
    "complement_data/{file}"
  shell:
    "cat {input}|tr atcg tagc > {output}"

rule reverse:
  input:
    "complement_data/{file}"
  output:
    "reversed_data/{file}"
  shell:
    "cat {input} |rev > {output}"
