#1.将json读进R里并简化
meta <- jsonlite::fromJSON("metadata.cart.2020-01-08.json")
colnames(meta)
#找了大半天的TCGA-ID就在associated_entities里，说出来你可能不信，data.frame的一个格子里竟然可以装列表，我今天也是第一次见。
entity <- meta$associated_entities
meta$associated_entities[[1]]
class(entity)

jh = function(x){
  as.character(x[3])
}
jh(entity[[1]])

class(jh(entity[[1]]))

xe = sapply(entity,jh)


#简化为file_name和TCGAid两列