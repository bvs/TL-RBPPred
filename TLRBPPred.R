TLRBPPred <- function(Xs,Xt,Ys,Yt,k,p,lambda,gama,iters) {
	print("Step: 1 Setting predefined variables: X, Y, m, c, ns, nt, YY");
	X <- as(cbind(Xs,Xt),"sparseMatrix") # data matrix of source and target domains
	Y <- as.matrix(rbind(as.matrix(Ys),as.matrix(Yt))) # label matrix of source and target domains
	m <- dim(X)[1] # number of features in shared feature space
	c <- length(unique(Y)) # number of classes in shared label space
	ns <- dim(Xs)[2] # number of samples in source
	nt <- dim(Xt)[2] # number of samples in target
	YY <- c()
	for(i in sort(as.vector(array(unique(as.vector(Y)), c(1, length(unique(Y))))))) {
		YY <- cbind(YY, as.vector(Y) == i)
		YY[YY == FALSE] <- 0
	}
	Y <- as.matrix(apply(YY,1,which.max))
	Ys <- YY[1:ns,];
	Yt <- YY[(ns+1):dim(YY)[1],];
	print("Step: 2 Normalization of Xs and Xt for Classification");
	Xs <-  Xs %*% as(diag(c(as.matrix(1/sqrt(apply(Xs^2,2,sum))))),"sparseMatrix")
	Xt <-  Xt %*% as(diag(c(as.matrix(1/sqrt(apply(Xt^2,2,sum))))),"sparseMatrix")
	print("Step: 3 Construction of Graph Laplacian for Xs and Xt (This will take a while!)");
	manifold.k <- p;
	manifold.Metric <- "Cosine";
	manifold.NeighborMode <- "KNN";
	manifold.WeightMode <- "Cosine";
	manifold.bNormalizeGraph <- 0;
	manifold.bNormalized <- 0;
	Wus <- affinity(as.matrix(Xs),manifold.k,manifold.Metric,manifold.NeighborMode,manifold.WeightMode,manifold.bNormalized);
	print("Wus created")
	Dus <- diagnal(Wus,manifold.bNormalizeGraph)
	Wut <- affinity(as.matrix(Xt),manifold.k,manifold.Metric,manifold.NeighborMode,manifold.WeightMode,manifold.bNormalized);
	print("Wut created")
	Dut <- diagnal(Wut,manifold.bNormalizeGraph)
	Wvt <- affinity(t(as.matrix(Xt)),manifold.k,manifold.Metric,manifold.NeighborMode,manifold.WeightMode,manifold.bNormalized);
	print("Wvt created")
	Dvt <- diagnal(Wvt,manifold.bNormalizeGraph)
	print("Step: 4 Initialization of variables: Us,Ut,H,Vs and Vt")
	Us <- matrix(runif(k*m),m) # k feature clusters in source domain
	Ut <- matrix(runif(k*m),m) # k feature clusters in target domain
	H <- t(matrix(runif(c*k),c)) # Association between U and V
	Vs <- 0.1+0.8*Ys; # c example classes in source domain
#	if(!is.null(Yt0) && length(Yt0$predictions) == nt) {
#		Vt <- c()
#		for(i in sort(as.vector(array(unique(as.vector(Yt0$predictions)), c(1, length(unique(Yt0$predictions))))))) {
#			Vt <- cbind(Vt, as.vector(Yt0$predictions) == i)
#			Vt[Vt == FALSE] <- 0
#		}
#		Yt0 <- c()
#		Vt <- 0.1+0.8*Vt; # c example classes in target domain 
#	} else {
		Vt <- matrix(runif(c*nt),nt)
#	}
	print(dim(Vt))
print(dim(H))
	print("Step: 5 Starting Graph Co-Regularized Nonnegative Matrix Tri-Factorization (GCMF) Algorithm");
	Acc <- c()
	Obj <- c()
	prob1 <- c()
	prob2 <- c()
	cls2 <- c()
	cls <- list()
	for(i in 1:iters) {
		Us <- Us * sqrt((Xs%*%(Vs%*%t(H))+lambda*Wus%*%Us)/(Us%*%(t(Us)%*%(Xs%*%Vs)%*%t(H))+lambda*Dus%*%Us+.Machine$double.eps))
		Us <- Us/((tcrossprod(rep(1,dim(Us)[1]),c(apply(Us^2,2,sum)^0.5)))+.Machine$double.eps)
	
		Ut <- Ut * sqrt((Xt%*%(Vt%*%t(H))+lambda*Wut%*%Ut)/(Ut%*%(t(Ut)%*%(Xt%*%Vt)%*%t(H))+lambda*Dut%*%Ut+.Machine$double.eps))
		Ut <- Ut/((tcrossprod(rep(1,dim(Ut)[1]),c(apply(Ut^2,2,sum)^0.5)))+.Machine$double.eps)
		
		Vs <- Vs/((tcrossprod(rep(1,dim(Vs)[1]),c(apply(Vs^2,2,sum)^0.5)))+.Machine$double.eps)

		Vt <- Vt * sqrt((t(Xt)%*%(Ut%*%H)+gama*Wvt%*%Vt)/(Vt%*%(t(Vt)%*%(t(Xt)%*%Ut)%*%H)+gama*Dvt%*%Vt+.Machine$double.eps))
		Vt <- Vt/((tcrossprod(rep(1,dim(Vt)[1]),c(apply(Vt^2,2,sum)^0.5)))+.Machine$double.eps)
		
		prob1 <- append(prob1, as.vector(Vt[,1]))
		prob2 <- append(prob2, as.vector(Vt[,2]))

		H <- H * sqrt((t(Us)%*%(Xs%*%Vs)+t(Ut)%*%(Xt%*%Vt))/(t(Us)%*%(Us%*%H%*%(t(Vs)%*%Vs))+t(Ut)%*%(Ut%*%H%*%(t(Vt)%*%Vt))+.Machine$double.eps))		

		Cls <- as.matrix(apply(Vt,1,which.max)) ### Predicted label based on highest probability in a column
		cls2 <- rbind(cls2,Cls)
		Lbl <- as.matrix(apply(Yt,1,which.max)) ### Ground truth label based on highest probability in a column
		acc <- (length(which(Cls == Lbl))/nt)*100; 
		Acc <- rbind(Acc,acc)
		O <- 0
		print(paste("iteration:",i, "Accuracy:",round(acc, digits=3)));
	}
	cls[['prob1']] <- prob1; cls[['prob2']] <- prob2; cls[['cls2']] <- cls2; cls[['Acc']] <- Acc;
	return(cls)
	print("Algorithm GCMF terminated!!!");
}
### End of GCMF ###
###
### affinity ###
affinity <- function(fea,k,Metric,NeighborMode,WeightMode,bNormalized) {
	library(Matrix)
	bSelfConnected <- 1;
	nSmp <- dim(fea)[1]
	maxM <- 62500000 # 500M
	BlockSize <- floor(maxM/(nSmp*3));
	if((NeighborMode == "KNN") && (k > 0)) {
		if(Metric == "Cosine") {
			if(!(bNormalized)) {
				#fea <- as(fea, "matrix");
				nSmp <- dim(fea)[1];
				nFea <- dim(fea)[2];
				if(!is.matrix(fea)) {
					fea2 <- t(fea);
					rm(fea);
					for(i in 1:nSmp) {
						fea2[,i] <- fea2[,i]/max(1e-10,sum(fea2[,i]^2)^0.5)
					}
					fea <- t(fea2);
					rm(fea2);
				} else { 
					feaNorm <- apply(fea^2,1,sum)^0.5;
					for(i in 1:nSmp) {
						fea[i,] <- fea[i,]/max(1e-12, feaNorm[i]);
					}
				}
			}
		}
		fea <- as(fea, "sparseMatrix");
		G <- array(0,c(nSmp*(k+1),3));
		for(i in 1:ceiling(nSmp/BlockSize)) {
			if(i == ceiling(nSmp/BlockSize)) {
				smpIdx <- ((i-1)*BlockSize+1):nSmp;
				dist <- fea[smpIdx,] %*% t(fea);
				dump <- t(apply(-dist,1,sort));
				idx <- t(apply(-dist,1,order));
				idx <- idx[,0:k+1];
				dump <- -dump[,0:k+1];
				G[((i-1)*BlockSize*(k+1)+1):(nSmp*(k+1)),1] <- kronecker(matrix(1,k+1,1),smpIdx);
				G[((i-1)*BlockSize*(k+1)+1):(nSmp*(k+1)),2] <- as.vector(idx);
				G[((i-1)*BlockSize*(k+1)+1):(nSmp*(k+1)),3] <- as.vector(dump);
			} else {
				smpIdx <- ((i-1)*BlockSize+1):(i*BlockSize);
				dist <- fea[smpIdx,] %*% t(fea);
				dump <- t(apply(-dist,1,sort));
				idx <- t(apply(-dist,1,order));
				idx <- idx[,0:k+1];
				dump <- -dump[,0:k+1];
				G[((i-1)*BlockSize*(k+1)+1):((i*BlockSize)*(k+1)),1] <- kronecker(matrix(1,k+1,1),smpIdx);
				G[((i-1)*BlockSize*(k+1)+1):((i*BlockSize)*(k+1)),2] <- as.vector(idx);
				G[((i-1)*BlockSize*(k+1)+1):((i*BlockSize)*(k+1)),3] <- as.vector(dump);
			}
		}
		library("Matrix")
#		W <- sparseMatrix(i=G[,1],j=G[,2],x=G[,3])
		W <- sparseMatrix(i=G[,1],j=G[,2],x=G[,3], dims=c(nSmp,nSmp))
		if(!(bSelfConnected)) {
       			for(i in 1:dim(W)[1]) {
				W[i,i] = 0;
			}
		}
		W <- pmax.sparse(W,t(W));
		return(W);
	}
}
### End of affinity ###
###
### pmax.sparse ### 
pmax.sparse <- function(..., na.rm = FALSE) {
	num.rows <- unique(sapply(list(...), nrow))
	num.cols <- unique(sapply(list(...), ncol))
	stopifnot(length(num.rows) == 1)
	stopifnot(length(num.cols) == 1)
	cat.summary <- do.call(rbind, lapply(list(...), summary))
	out.summary <- aggregate(x ~ i + j, data = cat.summary, max, na.rm)
	sparseMatrix(i = out.summary$i, j = out.summary$j, x = out.summary$x, dims = c(num.rows, num.cols))
}
### End of pmax.sparse ###
###
### diagnal ###
diagnal <- function(W, manifold.bNormalizeGraph) {
	if(manifold.bNormalizeGraph) {
		D <- as(1/sqrt(rowSums(W)), "spareMatrix");
		D[is.infinite(D)] <- 0;
		D <- diag(c(as.matrix(D)));
		W <- D %*% W %*% D;
		D <- as(D, "matrix");
		D[D>0] <- 1;
		D <- as(D,"sparseMatrix");
	} else {
		D <- as(diag(c(as.matrix(rowSums(W)))), "sparseMatrix");
		}
	return(D)
}
### End of diagnal ###
