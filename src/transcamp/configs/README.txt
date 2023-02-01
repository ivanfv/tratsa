Para meter la precisión de las distintas partes del algoritmo FlexFloat hay que
crear un archivo con el nombre de la serie y extensión .cfg, con cuatro líneas y
dos números por línea que serán el exponente y la mantisa de cada parte. La
serie quedará codificada con (covariance_exponent,covariance_mantissa).

################################################################################
covariance_exponent     covariance_mantissa
correlation_exponent    correlation_mantissa
stats_exponent          stats_mantissa
profile_exponent        profile_mantissa
################################################################################

