createOncoMatrix = function(maf){
     oncomat = data.table::dcast(data = maf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode, patientID)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                            fun.aggregate = function(x){
                              x = unique(as.character(x))
                              xad = x[x %in% c('Amp', 'Del')]
                              xvc = x[!x %in% c('Amp', 'Del')]

                              if(length(xvc)>0){
                                xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                              }

                              x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                              x = gsub(pattern = ';$', replacement = '', x = x)
                              x = gsub(pattern = '^;', replacement = '', x = x)
                              return(x)
                            } , value.var = 'Variant_Classification', fill = '')

    #If maf contains only one sample converting to matrix is not trivial.
    if(ncol(oncomat) == 2){
      genes = oncomat[,Hugo_Symbol]
      sampleId = colnames(oncomat)[2]
      oncomat = as.matrix(data.frame(row.names = genes, sample = oncomat[,2, with =FALSE]))
    }else if(nrow(oncomat) == 1){
      #If MAF has only one gene
      gene = oncomat[,Hugo_Symbol]
      oncomat[,Hugo_Symbol:= NULL]
      oncomat = as.matrix(oncomat)
      rownames(oncomat) = gene
      sampleID = colnames(oncomat)
      }else{
      oncomat = as.matrix(oncomat)
      rownames(oncomat) = oncomat[,1]
      oncomat = oncomat[,-1]
      }

     variant.classes = as.character(unique(maf[,Variant_Classification]))
     variant.classes = c('',variant.classes, 'Multi_Hit')
     names(variant.classes) = 0:(length(variant.classes)-1)

     #Complex variant classes will be assigned a single integer.
     vc.onc = unique(unlist(apply(oncomat, 2, unique)))
     vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
     names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
     variant.classes2 = c(variant.classes, vc.onc)

     oncomat.copy <- oncomat
    #Make a numeric coded matrix
    for(i in 1:length(variant.classes2)){
      oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
    }

    #If maf has only one gene
    if(nrow(oncomat) == 1){
      mdf  = t(matrix(as.numeric(oncomat)))
      rownames(mdf) = gene
      colnames(mdf) = sampleID
      return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
    }

    #convert from character to numeric
    mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
    rownames(mdf) = rownames(oncomat.copy)

    if(chatty){
      message('Sorting..')
    }

    #If MAF file contains a single sample, simple sorting is enuf.
    if(ncol(mdf) == 1){
      mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
      colnames(mdf) = sampleId

      oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
      colnames(oncomat.copy) = sampleId

      return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
    } else{
      #Sort by rows as well columns if >1 samples present in MAF
      #Add total variants per gene
      mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
        length(x[x != "0"])
      }))
      #Sort by total variants
      mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
      #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
      nMut = mdf[, ncol(mdf)]

      mdf = mdf[, -ncol(mdf)]

      mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

      mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
      tmdf = t(mdf) #transposematrix
      mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

      mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
      mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
      mdf = mdf.temp.copy

      #organise original character matrix into sorted matrix
      oncomat.copy <- oncomat.copy[,colnames(mdf)]
      oncomat.copy <- oncomat.copy[rownames(mdf),]

      return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
    }
}
