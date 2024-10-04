#-------------------------------------------------------------------------------------------------------------------------
split  = function(data, pi){
  data = data[sample(1:dim(data)[1], dim(data)[1]), ]
  index = (dim(data)[1]*pi)%>%as.integer()
  data_train = data[1:index,]
  data_test = data[(index + 1): dim(data)[1],]
  return(list(train = data_train, test = data_test))
}

#--------------------------------------------------------------------------------------------------------------------

Select_model_lasso = function(data_train){

  X = data_train[, -dim(data_train)[2]]
  y = data_train[, dim(data_train)[2]]

  # Séléction des variables les plus importantes
  lambda = cv.glmnet(as.matrix(X), as.matrix(y), alpha = 1)$lambda.1se

  model = glmnet(X, y , alpha = 1, lambda = lambda)

  vect_coef = as.matrix(coef(model)[-1]!= 0, nrow = 1)


  M_select = which(vect_coef!=0)

  return(list(m = M_select, vect_coef = vect_coef))
}


Select_model_aic = function(data_train, step){

  if (step  == "backward"){

    # Séléction des variables par AIC par backward elimination

    X = data_train[, -dim(data_train)[2]]
    y = data_train[, dim(data_train)[2]]
    data_train = data.frame(x= X, y = y)

    full_model = lm(y~., data = data_train)
    best_model = stepAIC(full_model, direction = "backward", k = log(dim(data_train)[1])) # k = log(n): = BIC
  }

  else if (step == "forward"){
    # Séléction des variables par AIC par backward elimination

    X = data_train[, -dim(data_train)[2]]
    y = data_train[, dim(data_train)[2]]
    data_train = data.frame(x = X, y = y)

    full_model = lm(y~., data = data_train)
    init_model = lm(y~ 1, data = data_train)
    best_model = stepAIC(init_model, scope = formula(full_model),  direction = "forward", k = log(dim(data_train)[1])) # k = log(n): = BIC
  }

  else{
    stop("Enter a valid direction")
  }

  # modèle selectionné
  M_init = colnames(data_train[, -dim(data_train)[2]])
  M_aic  = names(summary(best_model)$coefficients[, "Estimate"])[-1]

  #vect_coef = rep(0, length(M_init))
  #for (i in (1: length(M_init))){
  #  if (M_init[i]%in% M_aic){
  #    vect_coef[i] = 1
  #  }
  #}
  vect_coef = as.integer(M_init %in% M_aic)

  return( list(m = which(vect_coef!=0), vect_coef = vect_coef))
}

#-----------------------------------------------------------------------------------------------------------------------------------

# forcer une variable

matrices_MBV_lasso_parallel_forced = function (data, N, pi, numCores, indice){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1

  #index de la variable forcer

  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N)
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)
  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  # Mise en place de la parrallélisation
  cl = makePSOCKcluster(numCores)
  clusterExport(cl, varlist = c("pi", "p", "N", "indice",  "data", "split", "Select_model_lasso"), envir = environment())
  registerDoParallel(cl)

  mat = foreach(1:N, .combine = rbind, .multicombine = TRUE, .packages = c("glmnet", "magrittr")) %dopar%{

    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    s = Select_model_lasso(data_train)
    M_hat = s$m # les indexs du model selectionnés

    # forcer la variable
    M_hat = c(indice, M_hat)
    M_hat = M_hat%>%unique()


    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))

    vect_coef = s$vect_coef #stocke les coéfficients sélectionnés
    vect_coef[indice] = 1
    # récupère les estimations
    betaM_hat = model$coefficients[-1]
    intercept = model$coefficients[1]
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])^2



    return(list(vect_coef = vect_coef, betaM_hat = betaM_hat, intercept = intercept, index = M_hat))

  }
  stopCluster(cl)

  for (i in 1:N){
    matrice_M_hat[i, ] = unlist(mat[i, "vect_coef"])
    matriceBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "betaM_hat"])))
    #matriceVarBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "VarBetaM_hat"])))
  }
  vect_intercept = unname(unlist(mat[, "intercept"]))

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}


matrices_MBV_aic_parallel_forced = function (data, N, pi,  numCores, indice, direction){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1

  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N )
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)

  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  # Mise en place de la parrallélisation
  cl = makePSOCKcluster(numCores)
  clusterExport(cl, varlist = c("pi", "p", "N", "indice", "data", "split", "Select_model_aic"), envir = environment())
  registerDoParallel(cl)

  mat = foreach(1:N, .combine = rbind, .multicombine = TRUE, .packages = c("MASS", "magrittr")) %dopar%{

    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    select = Select_model_aic(data_train, step = direction)
    M_hat  = select$m

    #forcer la variable
    M_hat = c(indice, M_hat)
    M_hat = M_hat%>%unique()

    vect_coef = select$vect_coef
    vect_coef[indice] = 1

    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))
    # récupère les estimations
    betaM_hat = unname(model$coefficients[-1])
    intercept = unname(model$coefficients[1])
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])^2

    return(list(vect_coef = vect_coef, betaM_hat = betaM_hat, intercept = intercept, index = M_hat))

  }
  stopCluster(cl)

  for (i in 1:N){
    matrice_M_hat[i, ] = unlist(mat[i, "vect_coef"])
    matriceBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "betaM_hat"])))
    #matriceVarBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "VarBetaM_hat"])))
  }
  vect_intercept = unname(unlist(mat[, "intercept"]))

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}


matrices_MBV_lasso_parallel = function (data, N, pi, numCores){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1

  #index de la variable forcer

  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N)
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)
  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  # Mise en place de la parrallélisation
  cl = makePSOCKcluster(numCores)
  clusterExport(cl, varlist = c("pi", "p", "N", "data", "split", "Select_model_lasso"), envir = environment())
  registerDoParallel(cl)

  mat = foreach(1:N, .combine = rbind, .multicombine = TRUE, .packages = c("glmnet", "magrittr")) %dopar%{

    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    s = Select_model_lasso(data_train)
    M_hat = s$m # les indexs du model selectionnés


    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))

    vect_coef = s$vect_coef #stocke les coéfficients sélectionnés
    # récupère les estimations
    betaM_hat = model$coefficients[-1]
    intercept = model$coefficients[1]
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])^2



    return(list(vect_coef = vect_coef, betaM_hat = betaM_hat, intercept = intercept, index = M_hat))

  }
  stopCluster(cl)

  for (i in 1:N){
    matrice_M_hat[i, ] = unlist(mat[i, "vect_coef"])
    matriceBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "betaM_hat"])))
    #matriceVarBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "VarBetaM_hat"])))
  }
  vect_intercept = unname(unlist(mat[, "intercept"]))

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}

#-----------------------------------------------------------------------------------------

matrices_MBV_aic_parallel = function (data, N, pi,  numCores, direction){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1

  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N )
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)

  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  # Mise en place de la parrallélisation
  cl = makePSOCKcluster(numCores)
  clusterExport(cl, varlist = c("pi", "p", "N", "data", "split", "Select_model_aic"), envir = environment())
  registerDoParallel(cl)

  mat = foreach(1:N, .combine = rbind, .multicombine = TRUE, .packages = c("MASS", "magrittr")) %dopar%{

    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    select = Select_model_aic(data_train, step = direction)
    M_hat  = select$m
    vect_coef = select$vect_coef

    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))
    # récupère les estimations
    betaM_hat = unname(model$coefficients[-1])
    intercept = unname(model$coefficients[1])
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])^2

    return(list(vect_coef = vect_coef, betaM_hat = betaM_hat, intercept = intercept, index = M_hat))

  }
  stopCluster(cl)

  for (i in 1:N){
    matrice_M_hat[i, ] = unlist(mat[i, "vect_coef"])
    matriceBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "betaM_hat"])))
    #matriceVarBetaM_hat[i, unlist(mat[i, "index"])] = unname((unlist(mat[i, "VarBetaM_hat"])))
  }
  vect_intercept = unname(unlist(mat[, "intercept"]))

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}

#--------------------------------------------------------------------------------------------------------------------------------

matrices_MBV_lasso = function(data , N, pi){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1
  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N)
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)
  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  for (i in 1:N){
    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    s = Select_model_lasso(data_train)
    M_hat = s$m
    matrice_M_hat[i, ] = s$vect_coef #stocke les coéfficients sélectionnés

    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))

    # récupère les estimations
    betaM_hat = model$coefficients[-1]
    vect_intercept[i] = model$coefficients[1]
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])

    # stockage dans les matrices correspondantes
    matriceBetaM_hat[i, M_hat] = betaM_hat
    #matriceVarBetaM_hat[i, M_hat] = VarBetaM_hat
  }

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}

matrices_MBV_lasso_forced = function(data , N, pi, indice){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1
  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N)
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)
  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  for (i in 1:N){
    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    s = Select_model_lasso(data_train)
    M_hat = s$m
    M_hat = c(indice, M_hat)
    M_hat = M_hat%>%unique()

    #s$vect_coef [M_hat] = 1
    vect_coef = s$vect_coef
    vect_coef[indice] = 1
    matrice_M_hat[i, ] = vect_coef #stocke les coéfficients sélectionnés

    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))

    # récupère les estimations
    betaM_hat = model$coefficients[-1]
    vect_intercept[i] = model$coefficients[1]
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])

    # stockage dans les matrices correspondantes
    matriceBetaM_hat[i, M_hat] = betaM_hat
    #matriceVarBetaM_hat[i, M_hat] = VarBetaM_hat
  }

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}

matrices_MBV_aic = function(data , N, pi, direction){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1
  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N)
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)
  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  for (i in 1:N){
    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    s = Select_model_aic(data_train, step = direction)
    M_hat = s$m
    matrice_M_hat[i, ] = s$vect_coef #stocke les coéfficients sélectionnés

    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))

    # récupère les estimations
    betaM_hat = model$coefficients[-1]
    vect_intercept[i] = model$coefficients[1]
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])

    # stockage dans les matrices correspondantes
    matriceBetaM_hat[i, M_hat] = betaM_hat
    #matriceVarBetaM_hat[i, M_hat] = VarBetaM_hat
  }

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}

matrices_MBV_aic_forced = function(data , N, pi, indice, direction){

  # debuter l'enrigistrement
  tic.clearlog()
  tic()

  p = dim(data)[2] - 1
  # initialisation des différents matrices
  matrice_M_hat = matrix(0, ncol = p, nrow = N)
  matriceBetaM_hat = matrix(0, ncol = p, nrow = N)
  vect_intercept = rep(0, N)
  #matriceVarBetaM_hat = matrix(0, ncol = p, nrow = N)

  for (i in 1:N){
    s = split(data, pi)
    data_train = s$train
    data_test = s$test

    # selection du model sur le train
    s = Select_model_aic(data_train, step = direction)
    M_hat = s$m
    M_hat = c(indice, M_hat)
    M_hat = M_hat%>%unique()

    vect_coef = s$vect_coef
    vect_coef[indice] = 1
    matrice_M_hat[i, ] = vect_coef #stocke les coéfficients sélectionnés

    # data_select
    X = data_test[, M_hat]
    y = data_test[, dim(data_test)[2]]

    # estimations des coéfficients
    model = lm(y~., data = data.frame(X, y))

    # récupère les estimations
    betaM_hat = model$coefficients[-1]
    vect_intercept[i] = model$coefficients[1]
    #VarBetaM_hat = unname(summary(model)$coefficients[, "Std. Error"][-1])

    # stockage dans les matrices correspondantes
    matriceBetaM_hat[i, M_hat] = betaM_hat
    #matriceVarBetaM_hat[i, M_hat] = VarBetaM_hat
  }

  toc(log = TRUE, quiet = TRUE)
  logs = tic.log()
  time = as.numeric(gsub(" sec elapsed", "", logs[[1]]))
  tic.clearlog()

  return(list(m_M = matrice_M_hat, m_beta = matriceBetaM_hat, vect_intercept = vect_intercept, var_names = colnames(data[, -dim(data)[2]]), time = time))
}
#-----------------------------------------------------------------------------------------------------------------

vote2 = function(matrices, N, emp_alpha){

  p = dim(matrices$m_M)[2]

  duplications = matrices$m_M[duplicated(matrices$m_M),]

  unique = unique(duplications)

  compte = matrix(0, nrow = dim(unique)[1], ncol = 1) # stocke le nombre de vecteurs en duplications pour chaque vecteur dans unique
  for (i in 1:dim(unique)[1]){
    compte[i,] = sum(apply(duplications, MARGIN = 1, FUN = function(x) all(x == unique[i,])))
  }

  # récuperer l'index "ligne" de la matrice compte qui contient le plus de duplications

  index = which(compte == max(compte))[1] # on peut en avoir plusieurs

  M2_vote = unique[index,] == 1 # Transformer en TRUE et FALSE

  indice_k2 = which(apply(matrices$m_M , MARGIN = 1, FUN = function(x) all(x == M2_vote)))

  BM2_hat = matrices$m_beta[indice_k2, ]#M2_vote]
  #VarBeta_hat = matrices$m_Var[indice_k2, ]#M2_vote]

  intercept = matrices$vect_intercept[indice_k2]%>%mean() # intercept sur les lignes correspondant au modèle le plus fréquent

  # IC_empirique
  matrice_icEmp = data.frame(matrix(0, ncol = p, nrow = 2))
  colnames(matrice_icEmp) = paste0("x", 1:p)
  for (j in 1:p){
    matrice_icEmp[, j] = unname(quantile(BM2_hat[, j], c(emp_alpha/2, 1- (emp_alpha/2))))
  }


  beta2 = apply(BM2_hat, MARGIN = 2, FUN = mean)
  #Var2 = apply(VarBeta_hat, MARGIN = 2, FUN = mean)

  #IC_sup = qnorm((1+0.95)/2) * sqrt(Var2) + beta2
  #IC_inf = - qnorm((1+0.95)/2) * sqrt(Var2) + beta2
  ic_inf = unname(unlist(matrice_icEmp[1,]))
  ic_sup = unname(unlist(matrice_icEmp[2,]))
  visuel_ICemp  = data.frame(coef = matrices$var_names, ic_inf = ic_inf , ic_sup = ic_sup, beta = beta2, method = "ICemp")
  #visuel_IC = data.frame(coef = paste0("x",1:p), ic_inf = IC_inf, ic_sup = IC_sup, method = "IC")
  #data_visuel = rbind(visuel_IC, visuel_ICemp)

  gg = ggplot(data = visuel_ICemp, aes(x = coef))+
    geom_point(aes(y = beta, color = "beta"))+
    #geom_point(aes(y = beta_reel, color = "beta_real"))+
    geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup, color = "CI"), width = 0.2)+
    scale_color_manual(name = "", values = c("beta" = "red", "CI" = "darkgreen"))+
    labs(title = "Vote on the models", x = "coefs", y = "estimates")+
    theme_minimal()



  # visualisation de cette approche
  #gg = ggplot(data = data.frame(beta = beta2 , beta_reel = beta_j), aes(x =1:p))  +
  #geom_point(aes(y = beta, color = "beta_hat"))+
  #geom_point(aes(y = beta_reel, color = "beta_real"))+
  #geom_errorbar(data = data.frame(ic_inf = IC_inf, ic_sup = IC_sup), aes(ymin = ic_inf, ymax = ic_sup, color = "IC"), width = 0.2)+
  #scale_color_manual(name = "Legend", values = c("beta_real" = "black", "beta_hat" = "red", "IC" = "darkgreen"))+
  #labs(title = "vote on the most selected model",  x = "Coefs", y = "values")+
  #theme_minimal()
  ic = data.frame(matrice_icEmp)
  colnames(ic) = matrices$var_names

  return(list(gg = gg, m = M2_vote, beta = beta2, b0 = intercept, ic = ic, names_vars = matrices$var_names))
}

vote3= function(matrices, N, emp_alpha, s){

  p = dim(matrices$m_M)[2]

  M = matrices$m_M
  M_beta = matrices$m_beta
  #M_var = matrices$m_Var

  M_vote = apply(M, MARGIN = 2, FUN = mean) > s

  # recuperer les colonnes non nulle de M_vote
  indice_j = which(M_vote!=0)


  matrice_icEmp = matrix(0, ncol = p, nrow = 2)
  beta = rep(0, p)
  #var = rep(0, p)
  list = list() # initialiser un index_ligne qui va enrigistrer les index_ligne de chaque variable
  for (j in indice_j){
    #ligne non nulle
    index_ligne = which(M[,j]==1)
    list[[j]] = index_ligne
    beta[j] = mean(M_beta[index_ligne,j]) # récupèrer beta chapeau
    #var[j] = mean(M_var[index_ligne, j]) # recuperer les var chapeau
    matrice_icEmp[, j] = unname(quantile(M_beta[index_ligne, j], c(emp_alpha/2, 1- (emp_alpha/2)))) # IC empirique
  }

  # gerer l'intercept
  list = list[!sapply(list, is.null)] # enlever les index NULL
  names(list) = seq_along(list) # re-indexer la list
  vect_long = sapply(list, length)
  min_vect = list[[which.min(vect_long)]]

  intercept = matrices$vect_intercept[min_vect]%>%mean()

  # Comme decider de l'intercept ici???????



  # IC comme dans le vote1&2
  #IC_sup = qnorm((1+0.95)/2) * sqrt(var) + beta
  #IC_inf = - qnorm((1+0.95)/2) * sqrt(var) + beta

  M3_vote = beta!=0

  ic_inf = unname(unlist(matrice_icEmp[1,]))
  ic_sup = unname(unlist(matrice_icEmp[2,]))

  visuel_ICemp  = data.frame(coef = matrices$var_names, ic_inf = ic_inf , ic_sup = ic_sup, beta = beta , method = "ICemp")
  #visuel_IC = data.frame(coef = paste0("x",1:p), ic_inf = IC_inf, ic_sup = IC_sup, method = "IC")
  #data_visuel = rbind(visuel_IC, visuel_ICemp)

  gg = ggplot(data = visuel_ICemp , aes(x = coef))+
    geom_point(aes(y = beta, color = "beta"))+
    #geom_point(aes(y = beta_reel, color = "beta_real"))+
    geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup, color = "CI"), width = 0.2)+
    scale_color_manual(name = "", values = c("beta" = "red", "CI" = "darkgreen"))+
    labs(title = "Vote on the coefficients", x = "coefs", y = "estimates")+
    theme_minimal()


  #gg = ggplot( data = data.frame(beta_real = beta_j, beta_hat = beta ), aes(x = 1:p))+
  #geom_point(aes(y = beta_hat, color = "beta_hat"))+
  #geom_point(aes(y = beta_real, color = "beta_real"))+
  #geom_errorbar(data = data.frame(ic_inf = IC_inf, ic_sup = IC_sup), aes(ymin = ic_inf, ymax = ic_sup, color = "IC"), width = 0)+
  #geom_errorbar(data = data.frame(ic_inf = matrice_icEmp[1,], ic_sup = matrice_icEmp[2,]), aes(ymin = ic_inf, ymax = ic_sup, color = "IC empirique"), width = 0.2)+
  #scale_color_manual(name = "Legend", values = c("beta_real" = "black", "beta_hat" = "red", "IC empirique" ="blue", "IC"="darkgreen"))+
  #labs(title = "vote on the mean of columns coefs (ignoring lines where coef not selected)", x = "Coefs", y = "values")+
  #scale_y_continuous(breaks = seq(0, max(4.5), by = 0.5))+
  #theme_minimal()

  ic = data.frame(matrice_icEmp)
  colnames(ic) = matrices$var_names


  return(list(gg = gg, m = M3_vote, beta = beta , b0 = intercept, ic = ic, names_vars = matrices$var_names))
}


#' Function that handles storing our estimation and variable selection matrices during the different splits.
#'
#' @param method Method for variable selection. Should be one of \code{"Lasso"} or \code{"BIC"}.
#' @param cores Number of cores for parallel processing.
#' @param data Data set to be used for regression modeling.
#' @param formula Regression model to use, specified as a formula.
#' @param N Number of splits.
#' @param p_split Probabilities associated with the splits.
#' @param forced_var A character string specifying a predictor variable to be forced into selection. By default, it is NULL, allowing for no forced selection. If provided, this variable will be consistently selected during the N splits.
#' @param direction It can take two values: \code{"backward"} and \code{"forward"}. In the case of BIC, it specifies the direction in which the selection will be made.
#'
#' @details
#' We have data that we will split several times while shuffling it each time. Then, we will divide the data into two parts based on a specific probability for splitting. In the first half, we will perform model selection, followed by calibration on the second half. At the end of these steps, we will obtain matrices of dimensions N*p that represent the selected models and the estimated coefficients associated with these models.
#'
#' @importFrom stats coef lm quantile formula
#' @importFrom MASS stepAIC
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makePSOCKcluster clusterExport stopCluster
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar scale_color_manual labs theme_minimal
#' @importFrom magrittr `%>%`
#' @importFrom foreach foreach %dopar%
#' @import mlbench
#' @import tictoc
#' @return An object of class lmps
#'
#' @examples
#'
#' library(mlbench)
#' data("BostonHousing")
#' # lmps object
#' model = lmps(medv ~ ., data = BostonHousing, method = "Lasso", N = 50)
#'
#' \donttest{
#' # A parallelized example
#' # lmps object
#' model = lmps(medv ~ ., data = BostonHousing, method = "Lasso", N = 50, cores = 2)
#' }
#'
#' @export


lmps = function(formula, data, method, N, p_split = 0.5, cores = NULL, direction = "backward", forced_var = NULL) {
  # Vérification de la méthode
  if (!method %in% c("Lasso", "BIC")) {
    stop("Error on method argument, check ??lmps")
  }

  ##################  Regler la formule ##################################
  variables = all.vars(formula)
  target = variables[1]

  if ("." %in% all.names(formula)) {
    predictors = setdiff(names(data), target)
  } else {
    predictors = variables[-1]
  }

  data = data[, c(predictors, target), drop = FALSE]  # Assurez-vous que data reste un data.frame

  ######### Regler le problème forced_var ################
  if (!is.null(forced_var)) {
    if (!forced_var %in% predictors) {
      stop("Error: Enter a valid variable name for forced_var")
    }
    indice = which(predictors == forced_var)
  }

  ################# Regler la parallélisation ###################
  if (is.null(cores)) {
    if (method == "BIC") {
      matrices = if (is.null(forced_var)) {
        matrices_MBV_aic(data, N, pi = p_split, direction = direction)
      } else {
        matrices_MBV_aic_forced(data, N, pi = p_split, indice = indice, direction = direction)
      }
    } else {
      matrices = if (is.null(forced_var)) {
        matrices_MBV_lasso(data, N, pi = p_split)
      } else {
        matrices_MBV_lasso_forced(data, N, pi = p_split, indice = indice)
      }
    }
  } else {
    if (method == "BIC") {
      matrices = if (is.null(forced_var)) {
        matrices_MBV_aic_parallel(data, N, pi = p_split, numCores = cores, direction = direction)
      } else {
        matrices_MBV_aic_parallel_forced(data, N, pi = p_split, numCores = cores, indice = indice, direction = direction)
      }
    } else {
      matrices = if (is.null(forced_var)) {
        matrices_MBV_lasso_parallel(data, N, pi = p_split, numCores = cores)
      } else {
        matrices_MBV_lasso_parallel_forced(data, N, pi = p_split, indice = indice, numCores = cores)
      }
    }
  }

  # Créer la liste de résultats
  lmps = list(nb_splittings = N, matrix = matrices)
  class(lmps) = "lmps"

  return(invisible(lmps))
}


# lmps = function(formula, data, method , N, p_split = 0.5, cores = NULL, direction = "backward", forced_var = NULL){
#   # method = c("Lasso", "BIC")
#   # vote = c("model", "coef")
#   # numCores : = nombre de coeurs pour la parallélisation
#   # s.vect_coef := si jamais on utilise vote sur les coefficients, la proba de selection d'une variable
#   # data := data set
#   # formula : = model de regression
#   # N := nombres de splittages
#   # p_split := probabiltés de splittages
#
#   #### pour la méthode
#
#
#
#
#   if (method%in%c("Lasso", "BIC")==FALSE){
#     stop("Error on method argument, check ??lmps")
#   }
#
#
#
#   ##################  Regler la formule ##################################
#   variables = all.vars(formula)
#
#   target = variables[1]
#
#   if ("." %in% all.names(formula)) {
#     predictors = setdiff(names(data), target)
#   } else {
#
#     predictors = variables[-1]
#   }
#
#   data = data[, c(predictors, target)]
#
#   ######### Regler le problème forced_var
#
#   if(is.null(forced_var)){
#
#     ################# Regler la parrallélisation ###################
#     if(is.null(cores)){
#
#       # post_selection
#       if (method == "BIC"){
#         matrices = matrices_MBV_aic(data, N, pi = p_split, direction = direction)
#       }
#       else{
#         matrices = matrices_MBV_lasso(data, N, pi = p_split)
#       }
#
#     }
#
#     else{
#       # post_selection
#       if (method == "BIC"){
#         matrices = matrices_MBV_aic_parallel(data, N, pi = p_split, numCores = cores, direction = direction)
#       }
#       else{
#         matrices = matrices_MBV_lasso_parallel(data, N, pi = p_split, numCores = cores)
#       }
#
#
#     }
#
#   }
#
#   else{
#
#     if(forced_var %in% predictors){
#       indice = which(predictors == forced_var)
#
#       ################# Regler la parrallélisation ###################
#       if(is.null(cores)){
#
#         # post_selection
#         if (method == "BIC"){
#           matrices = matrices_MBV_aic_forced(data, N, pi = p_split, indice = indice, direction = direction)
#         }
#         else{
#           matrices = matrices_MBV_lasso_forced(data, N, pi = p_split, indice = indice)
#         }
#
#       }
#
#       else{
#         # post_selection
#         if (method == "BIC"){
#           matrices = matrices_MBV_aic_parallel_forced(data, N, pi = p_split, numCores = cores, indice = indice, direction = direction)
#         }
#         else{
#           matrices = matrices_MBV_lasso_parallel_forced(data, N, pi = p_split, indice = indice, numCores = cores)
#         }
#
#       }
#
#
#
#     }
#
#     else{
#       stop("Error: Entrer un nom de variable à forced valide")
#     }
#
#   }
#
#
#
#   lmps =  list(nb_splittings = N, matrix = matrices)
#   class(lmps) = "lmps"
#
#
#   return(invisible(lmps))
# }



# function summary

#' Summary function for our lmps object
#'
#'
#' @param object Our lmps object
#' @param ... Other arguments ignored (for compatibility with generic)
#' @return A summary of our lmps object
#'
#' @details
#' This function provides a summary of the data collected during the application of the lmps function.
#' It summarizes how many times the most frequently selected model was chosen across our N divisions, as well as
#' the selection frequency of variables in the different divisions.
#' It can also provide the execution time of the lmps function, which may vary significantly depending on the chosen post-selection method and the dimensionality of our data.
#'
#'
#'
#' @examples
#' library(mlbench)
#' data("BostonHousing")
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50)
#' summary(model)
#'
#'
#'\donttest{
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50, cores = 2)
#' summary(model)
#'}
#' @method summary lmps
#'@export
summary.lmps = function(object, ...){

  # temps d'éxécution
  cat("Execution time in seconds: ")
  cat("\n")
  print(object$matrix$time%>%unlist())
  cat("\n")
  cat("\n")

  # Nombre de slittages
  cat("Number of Splits")
  cat("\n")
  print(object$nb_splittings)
  cat("\n")
  cat("\n")
  # compter la fréquence des coefficients

  cat("Frequency of variable selection :")
  cat("\n")
  cat("\n")
  freq_coef = object$matrix$m_M%>%apply(MARGIN = 2, FUN = sum)%>%t()%>%data.frame()
  colnames(freq_coef) = object$matrix$var_names
  rownames(freq_coef) = c("Nb")

  freq_coef%>%print()

  cat("\n")
  cat("\n")
  # compter la repetition du model le plus frequent
  cat("The number of times the most frequent model was selected : ")
  cat("\n")

  matrices = object$matrix
  p = dim(matrices$m_M)[2]

  duplications = matrices$m_M[duplicated(matrices$m_M),]

  unique = unique(duplications)

  compte = matrix(0, nrow = dim(unique)[1], ncol = 1) # stocke le nombre de vecteurs en duplications pour chaque vecteur dans unique
  for (i in 1:dim(unique)[1]){
    compte[i,] = sum(apply(duplications, MARGIN = 1, FUN = function(x) all(x == unique[i,])))
  }

  compte%>%max()%>%print()

}


#-------------------------------------------------------------------------------

#' Creates an object of class CIps based on the provided parameters.
#'
#' @param x An object of class lmps, which contains the selection and coefficient estimation matrices.
#' @param vote The type of vote to perform: "model" for selection based on the most frequent model,
#' or "coef" for variable selection (e.g., if a variable is selected more than 50 percent of the time).
#' @param alpha Specifies the confidence level for the confidence intervals.
#' @param s.vote_coef A parameter between 0 and 1 that, when using "coef" voting,
#' indicates the frequency threshold for selecting a variable.
#' @return An object of class CIps.
#'
#' @details
#' After obtaining the lmps object, which provides the selection matrices (models and coefficients),
#' this function allows us to compute confidence intervals that are calculated empirically based on the chosen voting method and the desired level of certainty.
#' The confidence intervals are obtained through empirical calculation on each vector of estimates for the corresponding coefficients.
#'
#'CIps also provides an intercept (test version) estimated as follows: in the case of a vote on models, it takes the average of the intercept vector for the rows where the most frequently selected model in the N splits is chosen. For the vote on coefficients, the idea is to select the coefficient that has been chosen the least number of times among those retained and then average the intercept only for the rows where this coefficient is selected.
#'
#'@examples
#' library(mlbench)
#' data("BostonHousing")
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50)
#' # CIps object
#' cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
#'
#'
#'
#'\donttest{
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50, cores = 2)
#' # CIps object
#' cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
#'}
#' @export
#'
CIps = function(x, vote, alpha, s.vote_coef = 0.5){

  # vote error
  if (vote%in%c("model", "coef")==FALSE){
    stop("Error on vote argument, check ??CIps")
  }

  # maintenant qu'on est les matrices ==> passons au vote
  matrices = x$matrix
  N = x$nb_splittings

  if (vote == "model"){
    vote = vote2(matrices, N, emp_alpha = alpha)
  }
  else{
    vote = vote3(matrices, N, emp_alpha = alpha, s = s.vote_coef)
  }

  class(vote) = "CIps"

  return(invisible(vote))
}

## creer les methode print(), plot()

#' Print method for the CIps class
#'
#' It provides information on the selected variables, the estimated confidence intervals, and the coefficients of these selected variables.
#'
#' @param x An object of class CIps.
#' @param ... Additional arguments to be passed to the print function.
#'
#' @return No return value, called for its side effects, which is printing the object to the console.
#'
#'@examples
#' library(mlbench)
#' data("BostonHousing")
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50)
#' # CIps object
#' cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
#' # print
#' print(cips)
#'
#'\donttest{
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50, cores = 2)
#' # CIps object
#' cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
#' # print
#' print(cips)
#'}
#'
#'@method print CIps
#' @export
print.CIps = function(x, ...){
  cat("Variable selection : \n \n")
  select_vars = data.frame(x$m%>%t())
  colnames(select_vars) = c(x$names_vars)
  rownames(select_vars) = c("selected?")
  print(select_vars)
  cat("\n")
  cat("**************************************************************************************************************************************\n")
  cat('\n')
  cat("Estimates : \n")
  cat('\n')
  ic = x$ic
  #renvoyer que les intervalles de confiance des variables selectionnées
  ic_select = ic[, which(ic!=0, arr.ind = T)[,2]%>%unique]
  rownames(ic_select) = c("lower bound", "upper bound")
  beta = x$beta
  beta = data.frame(beta%>%t())
  colnames(beta) = c(x$names_vars)
  beta = beta[which(beta!=0)]
  rownames(beta) = c("beta")
  out = rbind(ic_select, beta)
  print(out)

}

#' Plot method for the CIps class
#'
#' It provides a ggplot graphic where the x-axis displays all the explanatory variables, with the confidence intervals of the selected variables shown in green and the coefficient estimates represented as red points.
#'
#' @param x An object of class CIps.
#' @param ... Additional arguments to be passed to the plot function.
#'
#' @return No return value, called for its side effects, which is plotting a graph.
#' @examples
#'
#' library(mlbench)
#' data("BostonHousing")
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50)
#' # CIps object
#' cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
#' # plot
#' plot(cips)
#'
#'\donttest{
#' # lmps object
#' model = lmps(medv~., data = BostonHousing, method = "Lasso", N = 50, cores = 2)
#' # CIps object
#'cips = CIps(model, vote = "coef", alpha = 0.05, s.vote_coef = 0.5)
#' # plot
#'plot(cips)
#'}
#'@method plot CIps
#' @export
plot.CIps = function(x, ...){
  # le ggplot
  plot(x$gg)
}


#' A predict function for cips
#'
#' This function generates predictions based on a cips object.
#'
#' @param object An object of class 'cips'.
#' @param newdata A dataframe containing new data to make predictions.
#' @param X Explanatory variables, default is NULL.
#' @param y Target corresponding to the explanatory variables, default is NULL.
#' @param ... Additional arguments for future use.
#'
#' @details When X and y are not NULL, they represent the data used to select the model.
#' These data are then reused by recalibrating on the selected subset to obtain the beta estimates (hybrid approach).
#' @return A numeric vector of predicted values.
#' @export
predict.CIps = function(object, newdata, X = NULL, y = NULL,  ...){


  if(dim(newdata)[2]!= length(object$names_vars)){
    stop("Error: different variable numbers")
  }

  # réarranger le data test comme names_vars
  X_test = newdata[, object$names_vars]
  X_test = X_test%>%as.matrix(ncol = length(object$numvars))%>%apply(MARGIN = 2, FUN = as.numeric)

  if(is.null(X)){

    pred = X_test%>%as.matrix()%*% as.matrix(object$beta) + object$b0
  }

  else{

    dim = dim(X)[2] #enrigistrer la dimension
    X = X[, object$names_vars] # réordonner au cas où
    vars_select = object$m%>%as.integer()
    X = X[, which(vars_select!=0)]
    data = data.frame(X = X, y = y)
    colnames(data) = c(paste0("x_", (1:dim(X)[2])), "y")
    model = lm(y~., data = data)
    b0 = model$coefficients[1]
    beta = rep(0, dim)
    beta [which(vars_select!=0)] = model$coefficients[-1]

    pred = X_test%>%as.matrix()%*%as.matrix(beta) + b0

  }


  return (pred)
}




