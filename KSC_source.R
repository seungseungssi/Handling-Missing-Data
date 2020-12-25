#####KSC 하는 함수 만들기###########

ksc <- function(data,k = 3, w= 0.9,max.iter = 10000){
    # divide the features into two parts
    F_o = data[,colSums(is.na(data))==0]
    F_m = data[,colSums(is.na(data))!=0]
    F_m[rowSums(is.na(F_m)) != 0,] <-NA
    
    # get distance matrix for s_i,j
    s_ij = dist(F_m, diag = T, upper = T)
    s_ij = as.matrix(s_ij) # s_ij is fixed no matter what
    s_ij[is.na(s_ij)] <- 0
    
    # initialize the cluster membership and centroids
    membership = sample(1:k, nrow(data),replace = T)
    initial_center = sample(1:nrow(data), k, replace = F)
    centroid = matrix(c(rep(0, k*ncol(F_o))), nrow = k)
    for (i in 1:k){
        centroid[i,] = F_o[initial_center[i],]
    }
    
    #centroid = matrix(c(rep(0, k*ncol(F_o))), nrow = k)
    #for (i in 1:k){
    #    centroid[i,]= apply(F_o[membership == i,], 2,mean)
    #}
    
    #get V_max for normalizing
    centroid_normal = apply(F_o,2,mean)
    V_max = sum(apply((F_o-centroid_normal)^2, 1,sum))
    
    #get CV_max for normalizing
    CV_max = sum(s_ij^2)
    
    # iter until convergence
    object_func = 0 # init
    for (i in 1:max.iter){
        ####get the arg min of cluster membership###
        #get the first term 
        temp1 = matrix(c(rep(0,k*nrow(data))), ncol = k )
        for (i in 1:nrow(F_o)){
            for (j in 1:k){
                center = centroid[j,]
                temp1[i,j] = sum((F_o[i,]-center)^2) 
            }
        }
        
        
        # get the second term 
        temp2 = matrix(c(rep(0,k*nrow(data))), ncol = k)
        for (i in 1:nrow(data)){
            backup = membership[i]
            for (j in 1:k){
                # i번쩨 obs가 j 번째 cluster로 배정될 겨우의 delta랑 second term을 구해준것임
                membership[i] = j
                delta = dist(membership, diag = T,upper=T)+10
                delta[delta == 10] <- 1
                delta[delta != 1] <- 0
                delta = as.matrix(delta)
                
                #get the second term
                delta_s_ij = 0 
                for (h in 1:nrow(data)){
                    delta_s_ij = delta_s_ij + delta[i,h]*s_ij[i,h]^2
                }
                temp2[i,j] = delta_s_ij 
            }
            membership[i] = which.min(w*temp1[i,]/V_max+(1-w)*temp2[i,]/CV_max)
        }
        
        # get the total objective function
        temp_sum = w*temp1/V_max + (1-w)*temp2/CV_max
        
        # update the membership
        #membership = apply(temp_sum,1, which.min)
        
        # get the first term in objective function
        new_object_func = sum(apply(temp_sum,1, min))
        
        # if converge break and return cluster membership 
        if (abs(new_object_func - object_func) <0.0001){
            return(list("membership" =  membership, "object_func"= new_object_func))
            break
        }
        
        #update the object_func
        else {object_func = new_object_func}
        
        #update the centroid
        for (i in 1:k){
            centroid[i,]= apply(F_o[membership == i,], 2, mean)
        }
        print(object_func)
    }
    print("did not converge")
    return(list("membership" =  membership, "object_func"= new_object_func))
}