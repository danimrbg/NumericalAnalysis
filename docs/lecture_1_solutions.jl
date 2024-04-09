# Questão 3a
expressao_questao_3a = absolute_error(Float64(pi), 22 / 7)
resposta_questao_3a = 0.0012644892673496777

# Questão 3b
expressao_questao_3b = absolute_error(Float64(pi), 3.1415)
resposta_questao_3b = 9.265358979293481e-5





# Questão 4a
expressao_questao_4a = relative_error(Float64(factorial(7)),sqrt(2*Float64(pi)*7)*(7/Float64(exp(1)))^7)
resposta_questao_4a = 0.011826223886416323

# Questão 4b
expressao_questao_4b = relative_error(Float64(factorial(20)),sqrt(2*Float64(pi)*20)*(20/Float64(exp(1)))^20)
resposta_questao_4b = 0.0041576526228797





# Questão 5a
expressao_questao_5a = approximated_value_interval_given_relative_error(Float64(10),1e-3)
resposta_questao_5a = (9.99, 10.01)

# Questão 5b
expressao_questao_5b = approximated_value_interval_given_relative_error(Float64(100),1e-3)
resposta_questao_5b = (99.9, 100.1)

# Questão 5c
expressao_questao_5c = approximated_value_interval_given_relative_error(Float64(-50),1e-3)
resposta_questao_5c = (-50.05, -49.95)

# Questão 5d
expressao_questao_5d = approximated_value_interval_given_relative_error(Float64(10),1e-4)
resposta_questao_5d = (9.999, 10.001)

# Questão 5e
expressao_questao_5e = approximated_value_interval_given_relative_error(Float64(100),1e-4)
resposta_questao_5e = (99.99, 100.01)

# Questão 5f
expressao_questao_5f = approximated_value_interval_given_relative_error(Float64(-50),1e-4)
resposta_questao_5f = (-50.005, -49.995)





# Questão 6
# Valor de x
x = (13 / 14 - 6 / 7) / (2 * exp(1) - 5.4)
valor_x = 1.9535401392860217
# Valor aproximado de x por truncamento na aritmética finita
trunc_x = trunc(trunc_sum(13 / 14, - 6 / 7, digits = 3) / trunc_sum(2 * exp(1), - 5.4, digits = 3), digits = 3)
valor_trunc_x = 2.028
# Valor aproximado de x por arredondamento na aritmética finita
round_x = round(round_sum(13 / 14, - 6 / 7, digits = 3) / round_sum(2 * exp(1), - 5.4, digits = 3), digits = 3)
valor_round_x = 1.946
# Erros absoluto e relativo quanto ao valor aproximado via truncamento
erro_abs_trunc = absolute_error(x, trunc_x)
valor_erro_abs_trunc = 0.07445986071397837
erro_rel_trunc = relative_error(x, trunc_x)
valor_erro_rel_trunc = 0.0381153472184052
# Erros absoluto e relativo quanto ao valor aproximado via arredondamento
erro_abs_round = absolute_error(x, round_x)
valor_erro_abs_round = 0.007540139286021708
erro_rel_round = relative_error(x, round_x)
valor_erro_rel_round = 0.0038597309235619145