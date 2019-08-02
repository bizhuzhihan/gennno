gennno
========

Gene annotation
--------

# CARD:
```
eg = gennno.dove(genelist, path)
eg.card(path_card)
eg.db
```
# 参数说明
genelist: 基因名称（list）
path: 工作路径，用于保存临时文件
path_card: card.json的路径，没有就不填，会自动下载到工作路径中
eg.db: 返回包含gene，description的pd.DataFrame

# UniProt
```
eg = gennno.dove(genelist, path)
eg.uniprot(species_name)
eg.db
```
# 参数说明
genelist: 基因名称（list）
path: 工作路径，用于保存临时文件
species_name: 物种名称，如 "Klebsiella pneumoniae"
eg.db: 返回包含gene，function的pd.DataFrame
