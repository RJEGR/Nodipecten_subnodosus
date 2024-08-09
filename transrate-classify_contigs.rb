# Proposito, calcular de manera independiente las metricas de componentes de alineamiento (classify.contigs)
# para resultados de multiples alineamiento
# esto debido a la obsolencia del uso de transrate en el server omica
def classify_contigs cutoff
  @assembly.classify_contigs cutoff
end

# Reduce all metrics for the assembly to a single quality score
# by taking the geometric mean of the scores for all contigs
# and multiplying it by the proportion of fragments whose most likely
# mapping is consistent with the assembly
# @return [Integer] the assembly score