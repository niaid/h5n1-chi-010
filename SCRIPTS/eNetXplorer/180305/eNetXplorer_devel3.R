library(glmnet)
library(calibrate)
require(gplots)
require(RColorBrewer)
library(progress)
#library(ROCR)

eNetXplorer <- function(x, ...) UseMethod("eNetXplorer")

# family = c("gaussian","binomial","poisson","multinomial","cox","mgaussian")???
eNetXplorer.default <- function(x, y, alpha=seq(0,1,by=0.2), family, nlambda=100, nlambda.ext=NULL,
seed=NULL, scaled=T, n_fold=5, n_run=25, n_perm_null=20, QF.FUN=NULL, QF_label="QF", cor_method="pearson", fold_distrib_fail.max=100, ...)
{
    if (family=="gaussian") {
        res <- eNetXplorerGaussian(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, nlambda.ext=nlambda.ext,seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,n_perm_null=n_perm_null, QF.FUN=QF.FUN, QF_label=QF_label, cor_method=cor_method, ...)
    } else if (family=="binomial") {
        res <- eNetXplorerBinomial(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, nlambda.ext=nlambda.ext, seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, QF.FUN=QF.FUN, QF_label=QF_label, fold_distrib_fail.max=fold_distrib_fail.max, ...)
    } else if (family=="multinomial") {
        res <- eNetXplorerMultinomial(x=x, y=y, alpha=alpha, family=family, nlambda=nlambda, nlambda.ext=nlambda.ext, seed=seed, scaled=scaled, n_fold=n_fold, n_run=n_run,
        n_perm_null=n_perm_null, QF.FUN=QF.FUN, QF_label=QF_label, fold_distrib_fail.max=fold_distrib_fail.max, ...)
    } else {
      stop("Family type not supported\n")
    }
    res$call <- match.call()
    class(res) <- "eNetXplorer"
    res
}

print.eNetXplorer <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    x$alpha_best_lambda
}

summary.eNetXplorer <- function(x, ...)
{
    coeff_mat = cbind(x$alpha,x$best_lambda,x$model_QF_est,x$QF_model_vs_null_pval)
    colnames(coeff_mat) = c("alpha","lambda.max","QF.est","model.vs.null.pval")
    res <- list(call=x$call,
    coefficients=coeff_mat
    )
    class(res) <- "summary.eNetXplorer"
    res
}

print.summary.eNetXplorer <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

pdf.summary <- function (x, path=NULL, filename=NULL, ...)
{
    if (is.null(path)) {
        path = getwd()
    }
    if (is.null(filename)) {
        filename = "eNetXplorer.summary.pdf"
    }
    pdf(file.path(path,filename),width=7,height=5)
    plot(x,plot.type="summary")
    for (i_alpha in 1:length(x$alpha)) {
        plot(x,plot.type="lambda.vs.QF",alpha.index=i_alpha)
        plot(x,plot.type="measured.vs.oob",alpha.index=i_alpha)
        if (x$family%in%c("binomial","multinomial")) {
            plot(x,plot.type="contingency",alpha.index=i_alpha)
        }
        for (stat in c("coeff","freq")) {
            plot(x,plot.type="feature.caterpillar",alpha.index=i_alpha,stat=stat)
            plot(x,plot.type="feature.heatmap",alpha.index=i_alpha,stat=stat)
        }
    }
    dev.off()
}

plot.eNetXplorer.summary <- function (
x, main=NULL, col.main="black", cex.main=0.95, show.pval.ref=T,
...)
{
    y1 = x$model_QF_est
    y2 = -log10(x$QF_model_vs_null_pval)
    
    par(mar = c(5,4,4,4)) # default is mar = c(lower=5,left=4,top=4,right=2) + 0.1.
    
    plot(x$alpha,y1,type="b",yaxt="n",col="red",xlab="alpha",ylab="")
    axis(side=2, las=3, cex.axis=0.85, col.axis="red")
    mtext(text="QF (response vs out-of-bag predicted)", side=2, col="red", line=2.1)
    
    if (max(y2)==min(y2)) {
        delta = 1
        baseline = min(y2)-0.5
    } else {
        delta = max(y2)-min(y2)
        baseline = min(y2)
    }
    rescale = (max(y1)-min(y1))/delta
    shift = min(y1)-baseline*rescale
    y2n = y2*rescale+shift
    
    lines(x$alpha,y2n,type="b",col="blue")
    axis(side=4, las=3, labels=signif((axTicks(2)-shift)/rescale,digits=3), at=axTicks(2),cex.axis=0.85, col.axis="blue")
    mtext(text="-log10 pval (model vs null)", side=4, col="blue", line=2.1)
    
    if (show.pval.ref==T) {
        pval_ref = c(0.001,0.01,0.05,0.1)
        pval_ref_transf = -log10(pval_ref)*rescale+shift
        par(mgp = c(0,-2.8,0)) # sets location of axis labels, tick mark labels, and tick marks.
        for (i in 1:length(pval_ref)) {
            abline(h=pval_ref_transf[i],lty=3,col="blue")
            axis(side=4, las=1, labels=paste0("pval=",pval_ref), at=pval_ref_transf+0.01*(max(y1)-min(y1)),cex.axis=0.6, col.axis="blue", tcl=0.) # tcl to set tick length and orientation
        }
        par(mgp=c(3,1,0)) # resets defaults
    }

    if (is.null(main)) {
        main = "model performance vs alpha values"
    }
    title(main, col.main=col.main, cex.main=cex.main)
}

plot.eNetXplorer.lambda.vs.QF <- function (
x, alpha_index_plot=alpha_index_plot, log="x", type="b",
xlab=NULL, ylab=NULL, main=NULL, col.main="black", cex.main=0.95,
...)
{
    if (is.null(xlab)) {
        xlab="lambda"
    }
    if (is.null(ylab)) {
        ylab="QF (response vs out-of-bag predicted)"
    }
    for (i_alpha in alpha_index_plot) {
        plot(x$lambda_values[[i_alpha]],x$lambda_QF_est[[i_alpha]],
        log=log,type=type,xlab=xlab,ylab=ylab,...)
        if (is.null(main)) {
            main = paste0("alpha=",x$alpha[i_alpha])
            if (x$QF_label!="QF") {
                main = paste0(main," ; QF=",x$QF_label)
            }
        }
        title(main, col.main = col.main, cex.main=cex.main)
        abline(v=x$best_lambda[i_alpha],lty=3)
    }
}

plot.eNetXplorer.contingency <- function (
x, alpha_index_plot=alpha_index_plot,
xlab=NULL, ylab=NULL, cex.lab=0.95, main=NULL, col.main = "black", cex.main=0.85,
cex.axis=1, symbol.size.inches=0.5, bg.color="steelblue2", fg.color=NULL, margin=0.2,
frequency.label=T, frequency.label.cex=1, frequency.label.offset=0,
...)
{
    if (is.null(xlab)) {
        xlab="class (true)"
    }
    if (is.null(ylab)) {
        ylab="class (predicted)"
    }
    for (i_alpha in alpha_index_plot) {
        class = levels(x$response)
        n_class = length(class)
        contingency = matrix(rep(NA,n_class**2),ncol=n_class) # rows=true, cols=predicted
        for (i_class in 1:n_class) {
            instance_select = x$response==class[i_class]
            contingency[i_class,] = apply(x$predicted_values[[i_alpha]][instance_select,],2,sum)
        }
        
        freq = NULL
        for (i in 1:nrow(contingency)) {
            for (j in 1:ncol(contingency)) {
                freq = rbind(freq,c(i,j,contingency[i,j]))
            }
        }
        
        x.min=1-(n_class-1)*margin
        x.max=n_class+(n_class-1)*margin
        symbols(x=freq[,1],y=freq[,2],circles=sqrt(freq[,3]/pi), inches=symbol.size.inches, bg=bg.color, fg=fg.color,xlab=xlab,ylab=ylab,xaxt="n",yaxt="n",panel.first = abline(h=1:n_class,v=1:n_class,lty=3),xlim=c(x.min,x.max),ylim=c(x.min,x.max))
        axis(1,at=1:n_class,labels=as.character(class),cex.axis=cex.axis)
        axis(2,at=1:n_class,labels=as.character(class),cex.axis=cex.axis)
        if (is.null(main)) {
            main = paste0("contingency matrix for alpha=",x$alpha[i_alpha])
        }
        title(main, col.main = col.main, cex.main=cex.main)
        if (frequency.label) {
            textxy(freq[,1],freq[,2],freq[,3],cex=frequency.label.cex,offset=frequency.label.offset)
        }
    }
}

plot.eNetXplorer.measured.vs.oob.categ <- function (
x, alpha_index_plot=alpha_index_plot,
xlab=NULL, ylab=NULL, cex.lab=0.95, main=NULL, col.main = "black", cex.main=0.85,
transparency=70, jitter=0.25, cex.pt=1.7, class.color=NULL,
instance.label=T, instance.label.cex = 0.45, instance.label.offset = 0.6,
...)
{
    n_instance = length(x$instance)
    if (is.null(xlab)) {
        xlab="response"
    }
    if (is.null(ylab)) {
        ylab="out-of-bag predicted (accuracy)"
    }
    for (i_alpha in alpha_index_plot) {
        accuracy = rep(NA,n_instance)
        for (i_instance in 1:n_instance) {
            accuracy[i_instance] = x$predicted_values[[i_alpha]][i_instance,as.numeric(x$response)[i_instance]]
        }
        
        boxplot(as.formula("accuracy ~ x.response"),data=data.frame(x$response,accuracy),
        xlab=xlab,ylab=ylab,cex.lab=cex.lab,outline=F,boxwex=0.5,range=0,...)
        if (is.null(main)) {
            main = paste0("alpha=",x$alpha[i_alpha],
            "; lambda=",signif(x$best_lambda[i_alpha],digits=3),
            "; QF=",signif(x$model_QF_est[i_alpha],digits=2))
        }
        title(main, col.main = col.main, cex.main=cex.main)
        
        class = levels(x$response)
        n_class = length(class)
        if (is.null(class.color)) {
            class.color = 1:n_class
        }
        for (i_class in 1:n_class) {
            instance_select = x$response==class[i_class]
            scatter = runif(sum(instance_select),-jitter,jitter)
            x_plot = i_class+scatter-mean(scatter)
            y_plot = accuracy[instance_select]
            col_par = as.numeric(col2rgb(class.color[i_class]))
            points(x_plot,y_plot,cex=cex.pt,col=rgb(col_par[1],col_par[2],col_par[3],
            transparency,maxColorValue=255),pch=16)
            if (instance.label) {
                textxy(x_plot,y_plot,x$instance[instance_select],
                cex=instance.label.cex,offset=instance.label.offset)
            }
        }
    }
}

plot.eNetXplorer.measured.vs.oob.numer <- function (
x, alpha_index_plot=alpha_index_plot,
instance.label=T,
instance.label.added.margin = 0.1,instance.label.cex = 0.5,instance.label.offset = 0.75,
xlab=NULL, ylab=NULL, cex.lab=0.95, main=NULL, col.main = "black", cex.main=0.85, col="red",
...)
{
    if (is.null(xlab)) {
        xlab="response"
    }
    if (is.null(ylab)) {
        ylab="out-of-bag predicted (mean, sd)"
    }
    for (i_alpha in alpha_index_plot) {
        
        pred_mean = x$predicted_values[[i_alpha]][,1]
        pred_sd = x$predicted_values[[i_alpha]][,2]
        
        x_range = range(x$response)
        if (!is.null(instance.label)&&is.logical(instance.label)&&instance.label==T) {
            x_range = x_range + c(-1,1)*(x_range[2]-x_range[1])*instance.label.added.margin
        }
        y_range = range(c(pred_mean-pred_sd,pred_mean+pred_sd),na.rm=T)
        
        plot(x_range,y_range,type="n",xlab=xlab,ylab=ylab,cex.lab=cex.lab,...)
        if (is.null(main)) {
            main = paste0("alpha=",x$alpha[i_alpha],
            "; lambda=",signif(x$best_lambda[i_alpha],digits=3),
            "; QF=",signif(x$model_QF_est[i_alpha],digits=2))
        }
        title(main, col.main = col.main, cex.main=cex.main)
        
        df = data.frame(x=x$response,y=pred_mean)
        mod <- lm(y ~ x, data = df)
        newx_range = x_range + c(-1,1)*(x_range[2]-x_range[1])*0.3
        newx <- seq(newx_range[1], newx_range[2], length.out=200)
        preds <- predict(mod, newdata = data.frame(x=newx),
        interval = 'confidence')
        polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'lightblue1', border = NA)
        abline(mod,lty=1,col="blue") # model
        lines(newx, preds[,3], lty=3, col="blue") # intervals
        lines(newx, preds[,2], lty=3, col="blue")
        
        points(x$response,pred_mean,pch=16,cex=1,col=col,lty=3,type="p")
        arrows(x$response,pred_mean-pred_sd,x$response,pred_mean+pred_sd,col=col,length=0.025,angle=90,code=3)
        if (instance.label) {
            textxy(x$response,pred_mean,x$instance,cex=instance.label.cex,offset=instance.label.offset)
        }
    }
}

plot.eNetXplorer.feature.heatmap <- function (
x, alpha_index_plot=alpha_index_plot, stat=stat,
feature.all = F,
feature.pval.thres = NULL, feature.set = NULL, feature.top.n = 25, signif.code = T,
col=NULL, breaks=NULL, xlab=NULL, ylab=NULL,
scale="none", Rowv=F, Colv=F, na.color=NULL,
cexRow = NULL, srtRow=0, cexCol=0.75, srtCol=45, margins=c(5,5),
key=T, dendogram="none", trace="none", main=NULL, col.main = "black", line = 1, cex.main=1.2,
notecol = "black", notecex = 1, subtitle1=NULL, col.subtitle1="black",line.subtitle1=-1,cex.subtitle1=0.65,
subtitle2=NULL, col.subtitle2="darkgray",line.subtitle2=-2,cex.subtitle2=0.55,
...) {
    n_alpha = length(x$alpha)
    for (i_alpha in alpha_index_plot) {
        if (stat=="freq") {
            feature_heatmap = as.matrix(x$feature_freq_mean)
            pval = x$feature_freq_model_vs_null_pval
            main.def = "Feature frequencies"
            if (x$family=="binomial") {
                class = levels(x$response)
            }
        } else if (stat=="coeff") {
            feature_heatmap = as.matrix(x$feature_coef_wmean)
            pval = x$feature_coef_model_vs_null_pval
            main.def = "Feature coefficients"
            if (x$family=="binomial") {
                class = levels(x$response)
                main.def = paste0(main.def," for class ",class[2])
            }
        }
        n_feature = length(x$feature)
        o = order(pval[,i_alpha])
        if (!feature.all) {
            if (!is.null(feature.pval.thres)) {
                n_feature = sum(pval[o,i_alpha]<feature.pval.thres)
            } else if (!is.null(feature.set)) {
                feature_select = which(x$feature%in%feature.set)
                n_feature = length(feature_select)
                o = feature_select[order(pval[feature_select,i_alpha])]
            } else {
                n_feature = min(n_feature,feature.top.n)
            }
        }
        if (n_feature>1) {
            o = o[1:n_feature]
            feature_heatmap = feature_heatmap[o,]
            colnames(feature_heatmap) = paste0("alpha = ",x$alpha)
            rownames(feature_heatmap) = x$feature[o]
            if (is.null(main)) {
                main = main.def
            }
            if (is.null(breaks)) {
                if (stat=="freq") {
                    breaks = seq(0,1,by=0.2)
                } else if (stat=="coeff") {
                    scale_max = quantile(abs(feature_heatmap),probs=0.95,na.rm=T)
                    scale_int = scale_max/5
                    scale_breaks = seq(scale_int,scale_max,by=scale_int)
                    breaks = c(-rev(scale_breaks),0,scale_breaks)
                }
            }
            if (is.null(col)) {
                if (stat=="freq") {
                    col = colorRampPalette(brewer.pal(length(breaks)-1,"Blues"))
                } else if (stat=="coeff") {
                    col = redgreen
                }
            }
            if (is.null(na.color)) {
                if (stat=="freq") {
                    na.color = "black"
                } else if (stat=="coeff") {
                    na.color = "white"
                }
            }
            if (is.null(cexRow)) {
                cexRow = min(1.6/log(nrow(feature_heatmap)),0.5)
            }
            
            if (signif.code) {
                x$feature_coef_model_vs_null_pval[o,]
                annot = matrix(rep("",n_feature*n_alpha),ncol=n_alpha)
                annot[which(pval[o,]<0.1,arr.ind=T)] = "."
                annot[which(pval[o,]<0.05,arr.ind=T)] = "*"
                annot[which(pval[o,]<0.01,arr.ind=T)] = "**"
                annot[which(pval[o,]<0.001,arr.ind=T)] = "***"
                heatmap.2(feature_heatmap, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color, col=col,
                          breaks = breaks,dendrogram = dendogram, margins=margins, cexRow=cexRow,
                          srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, key=key, trace=trace,
                          xlab= xlab, ylab = ylab,
                          cellnote = annot, notecol = notecol, notecex = notecex, ...)
                title(main, col.main=col.main, line=line, cex.main=cex.main)
                if (is.null(subtitle1)) {
                    subtitle1 = paste0("Features ranked by p-value for alpha=",x$alpha[i_alpha])
                }
                title(subtitle1,col.main=col.subtitle1, line=line.subtitle1, cex.main=cex.subtitle1)
                if (is.null(subtitle2)) {
                    subtitle2 = "P-value significance codes:  <0.001 (***), <0.01 (**), <0.05 (*), <0.1 (.)"
                }
                title(subtitle2,col.main=col.subtitle2, line=line.subtitle2, cex.main=cex.subtitle2)
            } else {
                heatmap.2(feature_heatmap, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color, col=col,
                breaks = breaks,dendrogram = dendogram, margins=margins, cexRow=cexRow,
                srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, key=key, trace=trace,
                xlab= xlab, ylab = ylab,
                ...)
                title(main, col.main=col.main, cex.main=cex.main)
            }
        }
    }
}

plot.eNetXplorer.feature.caterpillar <- function (
x, alpha_index_plot=alpha_index_plot, stat=stat, feature.all = F,
feature.pval.thres = NULL, feature.set = NULL, feature.top.n = 25, signif.code = T,
main = NULL, col.main = "black", line=1.5, cex.main=0.85,
subtitle=NULL, col.subtitle="darkgray", line.subtitle=0.5, cex.subtitle=0.55,
xlab=NULL, ylab=NULL, cex.lab=0.95,
...) {
    for (i_alpha in alpha_index_plot) {
        if (stat=="freq") {
            feature_mean = x$feature_freq_mean[,i_alpha]
            feature_sd = x$feature_freq_sd[,i_alpha]
            null_feature_mean = x$null_feature_freq_mean[,i_alpha]
            null_feature_sd = x$null_feature_freq_sd[,i_alpha]
            pval = x$feature_freq_model_vs_null_pval[,i_alpha]
            if (is.null(xlab)) {
                xlab = "feature frequency"
                if (x$family=="binomial") {
                    class = levels(x$response)
                }
            }
        } else if (stat=="coeff") {
            feature_mean = x$feature_coef_wmean[,i_alpha]
            feature_sd = x$feature_coef_wsd[,i_alpha]
            null_feature_mean = x$null_feature_coef_wmean[,i_alpha]
            null_feature_sd = x$null_feature_coef_wsd[,i_alpha]
            pval = x$feature_coef_model_vs_null_pval[,i_alpha]
            if (is.null(xlab)) {
                xlab = "feature coefficient"
                if (x$family=="binomial") {
                    class = levels(x$response)
                    xlab = paste0(xlab," for class ",class[2])
                }
            }
        }
        if (is.null(ylab)) {
            ylab = ""
        }
        if (is.null(main)) {
            if (stat=="freq") {
                plot_main = paste0("Feature frequencies ranked by p-value for alpha=",x$alpha[i_alpha])
            } else if (stat=="coeff") {
                plot_main = paste0("Feature coefficients ranked by p-value for alpha=",x$alpha[i_alpha])
            }
        } else {
            plot_main = main
        }
        n_feature = length(x$feature)
        o = order(pval)
        if (!feature.all) {
            if (!is.null(feature.pval.thres)) {
                n_feature = sum(pval[o]<feature.pval.thres)
            } else if (!is.null(feature.set)) {
                feature_select = which(x$feature%in%feature.set)
                n_feature = length(feature_select)
                o = feature_select[order(pval[feature_select])]
            } else {
                n_feature = min(n_feature,feature.top.n)
            }
        }
        if (n_feature>1) {
            o = rev(o[1:n_feature])
            feature = x$feature[o]
            feature_mean = feature_mean[o]
            feature_sd = feature_sd[o]
            null_feature_mean = null_feature_mean[o]
            null_feature_sd = null_feature_sd[o]
            pval = pval[o]
            x_range = range(c(feature_mean-feature_sd,feature_mean+feature_sd,
            null_feature_mean-null_feature_sd,null_feature_mean+null_feature_sd),na.rm=T)
            y_range = c(1,n_feature)
            if (stat=="freq") {
                x_range[1] = max(x_range[1],0)
                x_range[2] = min(x_range[2],1)
                if (x_range[1]==x_range[2]) { # alpha=0
                    x_range = x_range[1] + c(-1,1)*0.1
                }
            }
            plot(x_range,y_range,type="n",xlab=xlab,ylab=ylab,yaxt="n",cex.lab=cex.lab, ...)
            if (signif.code) {
                signif_legend = "P-value significance codes:  <0.001 (***), <0.01 (**), <0.05 (*), <0.1 (.)"
                signif_symbol = c("  .","  *"," **","***")
                signif_pval_thres = c(0.10,0.05,0.01,0.001)
                feature_signif = rep("   ",n_feature)
                for (i_signif in 1:length(signif_symbol)) {
                    feature_signif[pval<signif_pval_thres[i_signif]] = signif_symbol[i_signif]
                }
                feature = paste(feature,feature_signif)
            }
            axis(side=2,at=1:n_feature,labels=feature,las=2,cex.axis=min(1.6/log(n_feature),0.5),tck=-0.005)
            abline(h=(1:n_feature),col="lightgray",lty="dotted")
            if (signif.code) {
                if (is.null(subtitle)) {
                    subtitle = signif_legend
                }
                title(subtitle,col.main=col.subtitle, line=line.subtitle, cex.main=cex.subtitle)
            }
            title(plot_main, col.main=col.main, cex.main=cex.main)
            points(null_feature_mean,1:n_feature,pch=16,cex=1,col="blue",lty=3,type="p")
            arrows(pmax(x_range[1],null_feature_mean-null_feature_sd),1:n_feature,
            pmin(x_range[2],null_feature_mean+null_feature_sd),1:n_feature,col="blue",length=0.025,angle=90,code=3,lty=2)
            points(feature_mean,1:n_feature,pch=16,cex=1,col="red",lty=3,type="p")
            arrows(pmax(x_range[1],feature_mean-feature_sd),1:n_feature,
            pmin(x_range[2],feature_mean+feature_sd),1:n_feature,col="red",length=0.025,angle=90,code=3,lty=2)
            legend("topright",inset=c(0,-0.165),legend=c("model","null"),
            pch=16,col=c("red","blue"),xpd=T,cex=0.75)
        }
    }
}

plot.eNetXplorer <- function(x, plot.type, alpha.index=NULL, stat=NULL, ...)
#plot.type=c("lambda.vs.QF","measured.vs.oob","feature.heatmap","feature.caterpillar")
#stat=c("freq","coeff") if plot.type in c("feature.heatmap","feature.caterpillar")
{
    if (is.null(plot.type)|!plot.type%in%c("summary","lambda.vs.QF","measured.vs.oob","contingency",
    "feature.heatmap","feature.caterpillar")) {
        stop("parameter \'plot.type\' must be in c(\"summary\",\"lambda.vs.QF\",\"measured.vs.oob\"
        ,\"contingency\",\"feature.heatmap\",\"feature.caterpillar\")\n")
    }
    n_alpha = length(x$alpha)
    if (!is.null(alpha.index)) {
        if (is.numeric(alpha.index)&&all(floor(alpha.index)==alpha.index,na.rm=T) # check for integer
        &&alpha.index>0&&alpha.index<=n_alpha) {
            alpha_index_plot = alpha.index
        } else {
            stop("parameter \'alpha.index\' must be an integer in the range [1,length(alpha)]\n")
        }
    } else {
        alpha_index_plot = 1:n_alpha
    }
    
    if (plot.type=="summary") { # QF and model.vs.null.pval vs alpha
        plot.eNetXplorer.summary(x, ...)
    }
    
    if (plot.type=="lambda.vs.QF") { # QF distribution vs lambda
        plot.eNetXplorer.lambda.vs.QF(x, alpha_index_plot=alpha_index_plot, ...)
    }
    
    if (plot.type=="measured.vs.oob") { # true vs predicted response
        if (x$family=="gaussian") {
            plot.eNetXplorer.measured.vs.oob.numer(x, alpha_index_plot=alpha_index_plot, ...)
        } else if (x$family%in%c("binomial","multinomial")) {
            plot.eNetXplorer.measured.vs.oob.categ(x, alpha_index_plot=alpha_index_plot, ...)
        }
    }
    
    if (plot.type=="contingency") { # contingency matrix for categorical models
        if (x$family%in%c("binomial","multinomial")) {
            plot.eNetXplorer.contingency(x, alpha_index_plot=alpha_index_plot, ...)
        }
    }
    
    if (plot.type=="feature.heatmap") { # feature heatmaps
        if (!is.null(stat)&&stat%in%c("freq","coeff")) {
            plot.eNetXplorer.feature.heatmap (x, alpha_index_plot=alpha_index_plot, stat=stat, ...)
        } else {
            stop("parameter \'stat\' must be in c(\"freq\",\"coeff\")\n")
        }
    }
    
    if (plot.type=="feature.caterpillar") { # feature caterpillar plots
        if (!is.null(stat)&&stat%in%c("freq","coeff")) {
            plot.eNetXplorer.feature.caterpillar(x, alpha_index_plot=alpha_index_plot, stat=stat,...)
        } else {
            stop("parameter \'stat\' must be in c(\"freq\",\"coeff\")")
        }
    }
}

# Gaussian model
eNetXplorerGaussian <- function(x, y, alpha, family, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, QF.FUN, QF_label, cor_method, ...)
{
    n_instance = nrow(x)
    n_feature = ncol(x)
    instance = rownames(x)
    if (is.null(instance)) {
        instance = paste0("Inst.",1:nrow(x))
    }
    feature = colnames(x)
    if (is.null(feature)) {
        feature = paste0("Feat.",1:ncol(x))
    }
    prediction_type = "link"
    
    if (is.null(QF.FUN)) {
        QF = function(predicted,response) {
            cor_test = cor.test(predicted,response,method=cor_method)
            return (cor_test$estimate)
        }
        QF_label = paste0("correlation (",cor_method,")")
    } else {
        QF = QF.FUN
    }
    
    if (!is.null(seed)) {
        set.seed(seed) # set the seed for reproducibility
    }
    
    if (scaled) { # data are z-score transformed
        x_col_mean = apply(x,2,mean)
        x_col_sd = apply(x,2,sd)
        for (i in 1:ncol(x)) {
            x[,i] = (x[,i]-x_col_mean[i])/x_col_sd[i]
        }
    }
    
    foldid = NULL # we determine the fold template
    fold_size = floor(n_instance/n_fold)
    for (i_fold in 1:n_fold) {
        foldid = c(foldid,rep(i_fold,fold_size))
    }
    fold_rest = n_instance%%n_fold
    if (fold_rest>0) {
        for (i_fold in 1:fold_rest) {
            foldid = c(foldid,i_fold)
        }
    }
    
    n_alpha = length(alpha)
    best_lambda = rep(NA,n_alpha)
    model_QF_est = rep(NA,n_alpha)
    
    lambda_values = vector("list",n_alpha)
    lambda_QF_est = vector("list",n_alpha)
    
    foldid_per_run = vector("list",n_alpha)
    predicted_values = vector("list",n_alpha)
    
    feature_coef_wmean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_coef_wsd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_mean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_sd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    null_feature_coef_wmean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_coef_wsd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_freq_mean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_freq_sd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    QF_model_vs_null_pval = rep(NA,n_alpha)
    feature_coef_model_vs_null_pval = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_model_vs_null_pval = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    glmnet.control(factory=T) # resets internal glmnet parameters
    glmnet.control(...) # to allow changes of factory default parameters in glmnet
    
    for (i_alpha in 1:n_alpha) { # beginning of alpha loop
        fit = glmnet(x,y,alpha=alpha[i_alpha],family=family,nlambda=nlambda)
        # We extend the range of lambda values (if nlambda.ext is set)
        if (!is.null(nlambda.ext)&&nlambda.ext>nlambda) {
            lambda_max = fit$lambda[1]*sqrt(nlambda.ext/length(fit$lambda))
            lambda_min = fit$lambda[length(fit$lambda)]/sqrt(nlambda.ext/length(fit$lambda))
            lambda_values[[i_alpha]] = lambda_max*(lambda_min/lambda_max)**(((1:nlambda.ext)-1)/(nlambda.ext-1))
        } else {
            lambda_values[[i_alpha]] = fit$lambda
        }
        n_lambda = length(lambda_values[[i_alpha]])
        predicted_values_all_lambda = matrix(rep(NA,n_instance*n_run*n_lambda),ncol=n_lambda)
        feature_coef_per_lambda = vector("list",n_lambda)
        feature_freq_per_lambda = vector("list",n_lambda)
        for (i_lambda in 1:n_lambda) {
            feature_coef_per_lambda[[i_lambda]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
            feature_freq_per_lambda[[i_lambda]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
        }
        foldid_per_run[[i_alpha]] = matrix(rep(NA,n_instance*n_run),ncol=n_run)
        pb <- progress_bar$new(
        format = paste0("  running MODEL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            foldid_per_run[[i_alpha]][,i_run] = sample(foldid)
            feature_coef_per_run = vector("list",n_lambda)
            for (i_fold in 1:n_fold) {
                instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                fit = glmnet(x[instance_in_bag,],y[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]])
                predicted_values_all_lambda[n_instance*(i_run-1)+which(instance_out_of_bag),] = predict(fit, x[instance_out_of_bag,], type=prediction_type)
                for (i_lambda in 1:n_lambda) {
                    feature_coef_per_run[[i_lambda]] = cbind(feature_coef_per_run[[i_lambda]],coef(fit)[-1,][,i_lambda]) # we remove the intercept
                }
            }
            for (i_lambda in 1:n_lambda) {
                is.na(feature_coef_per_run[[i_lambda]]) <- feature_coef_per_run[[i_lambda]] == 0
                feature_coef_per_lambda[[i_lambda]][,i_run] = rowMeans(feature_coef_per_run[[i_lambda]], na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                feature_freq_per_lambda[[i_lambda]][,i_run] = rowSums(!is.na(feature_coef_per_run[[i_lambda]]))/n_fold
            }
            pb$tick()
        }
        
        lambda_QF_est_all_runs = matrix(rep(NA,n_lambda*n_run),ncol=n_run)
        for (i_lambda in 1:n_lambda) {
            for (i_run in 1:n_run) {
                lambda_QF_est_all_runs[i_lambda,i_run] = QF(predicted_values_all_lambda[(1:n_instance)+(i_run-1)*n_instance,i_lambda],y)
            }
        }
        lambda_QF_est[[i_alpha]] = apply(lambda_QF_est_all_runs,1,median)
        best_lambda_index = which.max(lambda_QF_est[[i_alpha]])
        
        best_lambda[i_alpha] = lambda_values[[i_alpha]][best_lambda_index]
        model_QF_est[i_alpha] = lambda_QF_est[[i_alpha]][best_lambda_index]
        
        predicted_values[[i_alpha]] = matrix(rep(NA,n_instance*2),ncol=2)
        for (i_instance in 1:n_instance) {
            predicted_values[[i_alpha]][i_instance,1] = median(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
            best_lambda_index])
            predicted_values[[i_alpha]][i_instance,2] = mad(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
            best_lambda_index])
        }
        
        weight_tot = rowSums(feature_freq_per_lambda[[best_lambda_index]])
        feature_coef_wmean[,i_alpha] = rowSums(feature_coef_per_lambda[[best_lambda_index]]*feature_freq_per_lambda[[best_lambda_index]],na.rm=T)/weight_tot
        for (i_feature in 1:n_feature) {
            feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(feature_freq_per_lambda[[best_lambda_index]][i_feature,]/weight_tot[i_feature] * (feature_coef_per_lambda[[best_lambda_index]][i_feature,] - feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
        }
        feature_freq_mean[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,mean)
        feature_freq_sd[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,sd)
        
        # permutation-based null
        null_QF_est_all_runs = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
        null_feature_coef = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        null_feature_freq = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        pb <- progress_bar$new(
        format = paste0("  running NULL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            for (i_perm_null in 1:n_perm_null) {   # permutation loop
                y_RDM = sample(y) # we randomly permute the response vector
                null_predicted_values = rep(NA,n_instance)
                null_feature_coef_per_run = NULL
                for (i_fold in 1:n_fold) {
                    instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                    instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                    fit = glmnet(x[instance_in_bag,],y_RDM[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]][best_lambda_index])
                    null_predicted_values[which(instance_out_of_bag)] = predict(fit, x[instance_out_of_bag,], s=lambda_values[[i_alpha]][best_lambda_index], type=prediction_type)
                    null_feature_coef_per_run = cbind(null_feature_coef_per_run,coef(fit)[-1,]) # we remove the intercept
                }
                null_QF_est_all_runs[i_run,i_perm_null] = QF(null_predicted_values,y_RDM)
                is.na(null_feature_coef_per_run) <- null_feature_coef_per_run == 0
                null_feature_coef[,(i_run-1)*n_perm_null+i_perm_null] = rowMeans(null_feature_coef_per_run, na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                null_feature_freq[,(i_run-1)*n_perm_null+i_perm_null] = rowSums(!is.na(null_feature_coef_per_run))/n_fold
            }
            pb$tick()
        }
        weight_tot = rowSums(null_feature_freq)
        null_feature_coef_wmean[,i_alpha] = rowSums(null_feature_coef*null_feature_freq,na.rm=T)/weight_tot
        for (i_feature in 1:n_feature) {
            null_feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(null_feature_freq[i_feature,]/weight_tot[i_feature] * (null_feature_coef[i_feature,] - null_feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
        }
        null_feature_freq_mean[,i_alpha] = apply(null_feature_freq,1,mean)
        null_feature_freq_sd[,i_alpha] = apply(null_feature_freq,1,sd)
        
        # Comparisons of model vs null
        n_tail = 0
        for (i_run in 1:n_run) {
            n_tail = n_tail + sum(null_QF_est_all_runs[i_run,]>lambda_QF_est_all_runs[best_lambda_index,i_run])
        }
        QF_model_vs_null_pval[i_alpha] = (n_tail+1)/(n_run*n_perm_null+1) # correction based on Phipson & Smyth (2010)
        
        for (i_feature in 1:n_feature) {
            n_coef = 0
            n_tail_coef = 0
            n_tail_freq = 0
            for (i_run in 1:n_run) {
                coef_null_vs_model = abs(null_feature_coef[i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)])>=abs(feature_coef_per_lambda[[best_lambda_index]][i_feature,i_run])
                n_coef = n_coef + sum(!is.na(coef_null_vs_model))
                n_tail_coef = n_tail_coef + sum(coef_null_vs_model,na.rm=T)
                n_tail_freq = n_tail_freq + sum(null_feature_freq[i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)]>=feature_freq_per_lambda[[best_lambda_index]][i_feature,i_run])
            }
            feature_coef_model_vs_null_pval[i_feature,i_alpha] = (n_tail_coef+1)/(n_coef+1)
            feature_freq_model_vs_null_pval[i_feature,i_alpha] = (n_tail_freq+1)/(n_run*n_perm_null+1)
        }
        
        cat("alpha = ",alpha[i_alpha]," completed.\n")
    } # end of alpha loop
    
    # return object
    list(
    # input data and parameters
    predictor = x, response = y, alpha = alpha, family = family, nlambda = nlambda,
    nlambda.ext = nlambda.ext, seed = seed, scaled = scaled, n_fold = n_fold, n_run = n_run,
    n_perm_null = n_perm_null, QF_label = QF_label, cor_method = cor_method, instance = instance,
    feature = feature, glmnet_params = glmnet.control(),
    # summary results
    best_lambda = best_lambda, model_QF_est = model_QF_est, QF_model_vs_null_pval = QF_model_vs_null_pval,
    # detailed results for plots and downstream analysis
    lambda_values = lambda_values, lambda_QF_est = lambda_QF_est,
    predicted_values = predicted_values,
    feature_coef_wmean = as(feature_coef_wmean,"CsparseMatrix"), feature_coef_wsd = as(feature_coef_wsd,"CsparseMatrix"),
    feature_freq_mean = as(feature_freq_mean,"CsparseMatrix"), feature_freq_sd = as(feature_freq_sd,"CsparseMatrix"),
    null_feature_coef_wmean = as(null_feature_coef_wmean,"CsparseMatrix"), null_feature_coef_wsd = as(null_feature_coef_wsd,"CsparseMatrix"),
    null_feature_freq_mean = as(null_feature_freq_mean,"CsparseMatrix"), null_feature_freq_sd = as(null_feature_freq_sd,"CsparseMatrix"),
    feature_coef_model_vs_null_pval = as(feature_coef_model_vs_null_pval,"CsparseMatrix"),
    feature_freq_model_vs_null_pval = as(feature_freq_model_vs_null_pval,"CsparseMatrix")
    
    )
}

# Binomial model
eNetXplorerBinomial <- function(x, y, alpha, family, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, QF.FUN, QF_label, fold_distrib_fail.max, ...)
{
    n_instance = nrow(x)
    n_feature = ncol(x)
    instance = rownames(x)
    feature = colnames(x)
    
    y <- as.factor(unlist(y)) # to force y into a factor
    class = levels(y)
    n_class = length(class)
    
    prediction_type = "class"
    
    if (is.null(QF.FUN)) {
        QF = function(predicted,response) {
            return (sum(predicted==response)/n_instance)
        }
        QF_label = "accuracy"
    }
    
    if (!is.null(QF.FUN)) {
        QF = QF.FUN
    }
    
    if (!is.null(seed)) {
        set.seed(seed) # set the seed for reproducibility
    }
    
    if (scaled) { # data are z-score transformed
        x_col_mean = apply(x,2,mean)
        x_col_sd = apply(x,2,sd)
        for (i in 1:ncol(x)) {
            x[,i] = (x[,i]-x_col_mean[i])/x_col_sd[i]
        }
    }
    
    foldid = NULL # we determine the fold template
    fold_size = floor(n_instance/n_fold)
    for (i_fold in 1:n_fold) {
        foldid = c(foldid,rep(i_fold,fold_size))
    }
    fold_rest = n_instance%%n_fold
    if (fold_rest>0) {
        for (i_fold in 1:fold_rest) {
            foldid = c(foldid,i_fold)
        }
    }
    
    n_alpha = length(alpha)
    best_lambda = rep(NA,n_alpha)
    model_QF_est = rep(NA,n_alpha)
    
    lambda_values = vector("list",n_alpha)
    lambda_QF_est = vector("list",n_alpha)
    
    foldid_per_run = vector("list",n_alpha)
    predicted_values = vector("list",n_alpha)
    
    feature_coef_wmean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_coef_wsd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_mean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_sd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    null_feature_coef_wmean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_coef_wsd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_freq_mean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_freq_sd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    QF_model_vs_null_pval = rep(NA,n_alpha)
    feature_coef_model_vs_null_pval = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_model_vs_null_pval = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    glmnet.control(factory=T) # resets internal glmnet parameters
    glmnet.control(...) # to allow changes of factory default parameters in glmnet
    
    for (i_alpha in 1:n_alpha) { # beginning of alpha loop
        fit = glmnet(x,y,alpha=alpha[i_alpha],family=family,nlambda=nlambda)
        # We extend the range of lambda values (if nlambda.ext is set)
        if (!is.null(nlambda.ext)&&nlambda.ext>nlambda) {
            lambda_max = fit$lambda[1]*sqrt(nlambda.ext/length(fit$lambda))
            lambda_min = fit$lambda[length(fit$lambda)]/sqrt(nlambda.ext/length(fit$lambda))
            lambda_values[[i_alpha]] = lambda_max*(lambda_min/lambda_max)**(((1:nlambda.ext)-1)/(nlambda.ext-1))
        } else {
            lambda_values[[i_alpha]] = fit$lambda
        }
        n_lambda = length(lambda_values[[i_alpha]])
        predicted_values_all_lambda = matrix(rep(NA,n_instance*n_run*n_lambda),ncol=n_lambda)
        feature_coef_per_lambda = vector("list",n_lambda)
        feature_freq_per_lambda = vector("list",n_lambda)
        for (i_lambda in 1:n_lambda) {
            feature_coef_per_lambda[[i_lambda]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
            feature_freq_per_lambda[[i_lambda]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
        }
        foldid_per_run[[i_alpha]] = matrix(rep(NA,n_instance*n_run),ncol=n_run)
        pb <- progress_bar$new(
        format = paste0("  running MODEL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            fold_distrib_OK = F
            fold_distrib_fail = 0
            while (!fold_distrib_OK) { # for categorical models, we must check that each class is in more than 1 fold.
                foldid_per_run[[i_alpha]][,i_run] = sample(foldid)
                if (sum(apply(table(foldid_per_run[[i_alpha]][,i_run],y)>0,2,sum)>1)==n_class) {
                    fold_distrib_OK = T
                } else {
                    fold_distrib_fail = fold_distrib_fail + 1
                    if (fold_distrib_fail>fold_distrib_fail.max) {
                        stop("Failure in class distribution across cross-validation folds. Increase n_fold, remove very small classes, or re-run with larger fold_distrib_fail.max\n")
                    }
                }
            }
            feature_coef_per_run = vector("list",n_lambda)
            for (i_fold in 1:n_fold) {
                instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                fit = glmnet(x[instance_in_bag,],y[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]])
                predicted_values_all_lambda[n_instance*(i_run-1)+which(instance_out_of_bag),] = predict(fit, x[instance_out_of_bag,], type=prediction_type)
                for (i_lambda in 1:n_lambda) {
                    feature_coef_per_run[[i_lambda]] = cbind(feature_coef_per_run[[i_lambda]],coef(fit)[-1,][,i_lambda]) # we remove the intercept
                }
            }
            for (i_lambda in 1:n_lambda) {
                is.na(feature_coef_per_run[[i_lambda]]) <- feature_coef_per_run[[i_lambda]] == 0
                feature_coef_per_lambda[[i_lambda]][,i_run] = rowMeans(feature_coef_per_run[[i_lambda]], na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                feature_freq_per_lambda[[i_lambda]][,i_run] = rowSums(!is.na(feature_coef_per_run[[i_lambda]]))/n_fold
            }
            pb$tick()
        }
        
        lambda_QF_est_all_runs = matrix(rep(NA,n_lambda*n_run),ncol=n_run)
        for (i_lambda in 1:n_lambda) {
            for (i_run in 1:n_run) {
                lambda_QF_est_all_runs[i_lambda,i_run] = QF(predicted_values_all_lambda[(1:n_instance)+(i_run-1)*n_instance,i_lambda],y)
            }
        }
        lambda_QF_est[[i_alpha]] = apply(lambda_QF_est_all_runs,1,median)
        best_lambda_index = which.max(lambda_QF_est[[i_alpha]])
        
        best_lambda[i_alpha] = lambda_values[[i_alpha]][best_lambda_index]
        model_QF_est[i_alpha] = lambda_QF_est[[i_alpha]][best_lambda_index]
        
        predicted_values[[i_alpha]] = matrix(rep(NA,n_instance*n_class),ncol=n_class)
        for (i_instance in 1:n_instance) {
            for (i_class in 1:n_class) {
                predicted_values[[i_alpha]][i_instance,i_class] =  sum(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
                best_lambda_index]==class[i_class])/n_run
            }
        }
        
        weight_tot = rowSums(feature_freq_per_lambda[[best_lambda_index]])
        feature_coef_wmean[,i_alpha] = rowSums(feature_coef_per_lambda[[best_lambda_index]]*feature_freq_per_lambda[[best_lambda_index]],na.rm=T)/weight_tot
        for (i_feature in 1:n_feature) {
            feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(feature_freq_per_lambda[[best_lambda_index]][i_feature,]/weight_tot[i_feature] * (feature_coef_per_lambda[[best_lambda_index]][i_feature,] - feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
        }
        feature_freq_mean[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,mean)
        feature_freq_sd[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,sd)
        
        # permutation-based null
        null_QF_est_all_runs = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
        null_feature_coef = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        null_feature_freq = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        pb <- progress_bar$new(
        format = paste0("  running NULL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            for (i_perm_null in 1:n_perm_null) {   # permutation loop
                fold_distrib_OK = F
                fold_distrib_fail = 0
                while (!fold_distrib_OK) { # for categorical models, we must check that each class is in more than 1 fold.
                    y_RDM = sample(y) # we randomly permute the response vector
                    if (sum(apply(table(foldid_per_run[[i_alpha]][,i_run],y_RDM)>0,2,sum)>1)==n_class) {
                        fold_distrib_OK = T
                    } else {
                        fold_distrib_fail = fold_distrib_fail + 1
                        if (fold_distrib_fail>fold_distrib_fail.max) {
                            stop("Failure in class distribution across cross-validation folds. Increase n_fold, remove very small classes, or re-run with larger fold_distrib_fail.max\n")
                        }
                    }
                }
                null_predicted_values = rep(NA,n_instance)
                null_feature_coef_per_run = NULL
                for (i_fold in 1:n_fold) {
                    instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                    instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                    fit = glmnet(x[instance_in_bag,],y_RDM[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]][best_lambda_index])
                    null_predicted_values[which(instance_out_of_bag)] = predict(fit, x[instance_out_of_bag,], s=lambda_values[[i_alpha]][best_lambda_index], type=prediction_type)
                    null_feature_coef_per_run = cbind(null_feature_coef_per_run,coef(fit)[-1,]) # we remove the intercept
                }
                null_QF_est_all_runs[i_run,i_perm_null] = QF(null_predicted_values,y_RDM)
                is.na(null_feature_coef_per_run) <- null_feature_coef_per_run == 0
                null_feature_coef[,(i_run-1)*n_perm_null+i_perm_null] = rowMeans(null_feature_coef_per_run, na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                null_feature_freq[,(i_run-1)*n_perm_null+i_perm_null] = rowSums(!is.na(null_feature_coef_per_run))/n_fold
            }
            pb$tick()
        }
        weight_tot = rowSums(null_feature_freq)
        null_feature_coef_wmean[,i_alpha] = rowSums(null_feature_coef*null_feature_freq,na.rm=T)/weight_tot
        for (i_feature in 1:n_feature) {
            null_feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(null_feature_freq[i_feature,]/weight_tot[i_feature] * (null_feature_coef[i_feature,] - null_feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
        }
        null_feature_freq_mean[,i_alpha] = apply(null_feature_freq,1,mean)
        null_feature_freq_sd[,i_alpha] = apply(null_feature_freq,1,sd)
        
        # Comparisons of model vs null
        n_tail = 0
        for (i_run in 1:n_run) {
            n_tail = n_tail + sum(null_QF_est_all_runs[i_run,]>lambda_QF_est_all_runs[best_lambda_index,i_run])
        }
        QF_model_vs_null_pval[i_alpha] = (n_tail+1)/(n_run*n_perm_null+1) # correction based on Phipson & Smyth (2010)
        
        for (i_feature in 1:n_feature) {
            n_coef = 0
            n_tail_coef = 0
            n_tail_freq = 0
            for (i_run in 1:n_run) {
                coef_null_vs_model = abs(null_feature_coef[i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)])>=abs(feature_coef_per_lambda[[best_lambda_index]][i_feature,i_run])
                n_coef = n_coef + sum(!is.na(coef_null_vs_model))
                n_tail_coef = n_tail_coef + sum(coef_null_vs_model,na.rm=T)
                n_tail_freq = n_tail_freq + sum(null_feature_freq[i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)]>=feature_freq_per_lambda[[best_lambda_index]][i_feature,i_run])
            }
            feature_coef_model_vs_null_pval[i_feature,i_alpha] = (n_tail_coef+1)/(n_coef+1)
            feature_freq_model_vs_null_pval[i_feature,i_alpha] = (n_tail_freq+1)/(n_run*n_perm_null+1)
        }
        
        cat("alpha = ",alpha[i_alpha]," completed.\n")
    } # end of alpha loop
    
    # return object
    list(
    # input data and parameters
    predictor = x, response = y, alpha = alpha, family = family, nlambda = nlambda,
    nlambda.ext = nlambda.ext, seed = seed, scaled = scaled, n_fold = n_fold, n_run = n_run,
    n_perm_null = n_perm_null, QF_label = QF_label, instance = instance,
    feature = feature, fold_distrib_fail.max = fold_distrib_fail.max, glmnet_params = glmnet.control(),
    # summary results
    best_lambda = best_lambda, model_QF_est = model_QF_est, QF_model_vs_null_pval = QF_model_vs_null_pval,
    # detailed results for plots and downstream analysis
    lambda_values = lambda_values, lambda_QF_est = lambda_QF_est,
    predicted_values = predicted_values,
    feature_coef_wmean = as(feature_coef_wmean,"CsparseMatrix"), feature_coef_wsd = as(feature_coef_wsd,"CsparseMatrix"),
    feature_freq_mean = as(feature_freq_mean,"CsparseMatrix"), feature_freq_sd = as(feature_freq_sd,"CsparseMatrix"),
    null_feature_coef_wmean = as(null_feature_coef_wmean,"CsparseMatrix"), null_feature_coef_wsd = as(null_feature_coef_wsd,"CsparseMatrix"),
    null_feature_freq_mean = as(null_feature_freq_mean,"CsparseMatrix"), null_feature_freq_sd = as(null_feature_freq_sd,"CsparseMatrix"),
    feature_coef_model_vs_null_pval = as(feature_coef_model_vs_null_pval,"CsparseMatrix"),
    feature_freq_model_vs_null_pval = as(feature_freq_model_vs_null_pval,"CsparseMatrix")
    )
}

# Multinomial model
eNetXplorerMultinomial <- function(x, y, alpha, family, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, QF.FUN, QF_label, fold_distrib_fail.max,...)
{
    n_instance = nrow(x)
    n_feature = ncol(x)
    instance = rownames(x)
    feature = colnames(x)
    
    y <- as.factor(unlist(y)) # to force y into a factor
    class = levels(y)
    n_class = length(class)
    
    prediction_type = "class"
    
    if (is.null(QF.FUN)) {
        QF = function(predicted,response) {
            avg_acc = 0
            for (i_class in 1:n_class) {
                avg_acc = avg_acc + sum((response==class[i_class])==(predicted==class[i_class]))
            }
            return (avg_acc/(n_class*n_instance))
        }
        QF_label = "average accuracy"
    } else {
        QF = QF.FUN
    }
    
    if (!is.null(seed)) {
        set.seed(seed) # set the seed for reproducibility
    }
    
    if (scaled) { # data are z-score transformed
        x_col_mean = apply(x,2,mean)
        x_col_sd = apply(x,2,sd)
        for (i in 1:ncol(x)) {
            x[,i] = (x[,i]-x_col_mean[i])/x_col_sd[i]
        }
    }
    
    foldid = NULL # we determine the fold template
    fold_size = floor(n_instance/n_fold)
    for (i_fold in 1:n_fold) {
        foldid = c(foldid,rep(i_fold,fold_size))
    }
    fold_rest = n_instance%%n_fold
    if (fold_rest>0) {
        for (i_fold in 1:fold_rest) {
            foldid = c(foldid,i_fold)
        }
    }
    
    n_alpha = length(alpha)
    best_lambda = rep(NA,n_alpha)
    model_QF_est = rep(NA,n_alpha)
    
    lambda_values = vector("list",n_alpha)
    lambda_QF_est = vector("list",n_alpha)
    
    foldid_per_run = vector("list",n_alpha)
    predicted_values = vector("list",n_alpha)
    
    feature_coef_wmean = vector("list",n_class)
    feature_coef_wsd = vector("list",n_class)
    feature_freq_mean = vector("list",n_class)
    feature_freq_sd = vector("list",n_class)
    for (i_class in 1:n_class) {
        feature_coef_wmean[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        feature_coef_wsd[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        feature_freq_mean[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        feature_freq_sd[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    }
    null_feature_coef_wmean = vector("list",n_class)
    null_feature_coef_wsd = vector("list",n_class)
    null_feature_freq_mean = vector("list",n_class)
    null_feature_freq_sd = vector("list",n_class)
    for (i_class in 1:n_class) {
        null_feature_coef_wmean[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        null_feature_coef_wsd[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        null_feature_freq_mean[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        null_feature_freq_sd[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    }
    
    QF_model_vs_null_pval = rep(NA,n_alpha)
    feature_coef_model_vs_null_pval =  vector("list",n_class)
    feature_freq_model_vs_null_pval =  vector("list",n_class)
    for (i_class in 1:n_class) {
        feature_coef_model_vs_null_pval[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
        feature_freq_model_vs_null_pval[[i_class]] = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    }
    
    glmnet.control(factory=T) # resets internal glmnet parameters
    glmnet.control(...) # to allow changes of factory default parameters in glmnet
    
    for (i_alpha in 1:n_alpha) { # beginning of alpha loop
        fit = glmnet(x,y,alpha=alpha[i_alpha],family=family,nlambda=nlambda)
        # We extend the range of lambda values (if nlambda.ext is set)
        if (!is.null(nlambda.ext)&&nlambda.ext>nlambda) {
            lambda_max = fit$lambda[1]*sqrt(nlambda.ext/length(fit$lambda))
            lambda_min = fit$lambda[length(fit$lambda)]/sqrt(nlambda.ext/length(fit$lambda))
            lambda_values[[i_alpha]] = lambda_max*(lambda_min/lambda_max)**(((1:nlambda.ext)-1)/(nlambda.ext-1))
        } else {
            lambda_values[[i_alpha]] = fit$lambda
        }
        n_lambda = length(lambda_values[[i_alpha]])
        predicted_values_all_lambda = matrix(rep(NA,n_instance*n_run*n_lambda),ncol=n_lambda)
        feature_coef_per_lambda = vector("list",n_lambda)
        feature_freq_per_lambda = vector("list",n_lambda)
        for (i_lambda in 1:n_lambda) {
            feature_coef_per_lambda[[i_lambda]] = vector("list",n_class)
            feature_freq_per_lambda[[i_lambda]] = vector("list",n_class)
            for (i_class in 1:n_class) {
                feature_coef_per_lambda[[i_lambda]][[i_class]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
                feature_freq_per_lambda[[i_lambda]][[i_class]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
            }
        }
        foldid_per_run[[i_alpha]] = matrix(rep(NA,n_instance*n_run),ncol=n_run)
        pb <- progress_bar$new(
        format = paste0("  running MODEL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            fold_distrib_OK = F
            fold_distrib_fail = 0
            while (!fold_distrib_OK) { # for categorical models, we must check that each class is in more than 1 fold.
                foldid_per_run[[i_alpha]][,i_run] = sample(foldid)
                if (sum(apply(table(foldid_per_run[[i_alpha]][,i_run],y)>0,2,sum)>1)==n_class) {
                    fold_distrib_OK = T
                } else {
                    fold_distrib_fail = fold_distrib_fail + 1
                    if (fold_distrib_fail>fold_distrib_fail.max) {
                        stop("Failure in class distribution across cross-validation folds. Increase n_fold, remove very small classes, or re-run with larger fold_distrib_fail.max\n")
                    }
                }
            }
            feature_coef_per_run = vector("list",n_lambda)
            for (i_lambda in 1:n_lambda) {
                feature_coef_per_run[[i_lambda]] = vector("list",n_class)
            }
            for (i_fold in 1:n_fold) {
                instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                fit = glmnet(x[instance_in_bag,],y[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]])
                predicted_values_all_lambda[n_instance*(i_run-1)+which(instance_out_of_bag),] = predict(fit, x[instance_out_of_bag,], type=prediction_type)
                for (i_lambda in 1:n_lambda) {
                    for (i_class in 1:n_class) {
                        feature_coef_per_run[[i_lambda]][[i_class]] = cbind(feature_coef_per_run[[i_lambda]][[i_class]],coef(fit)[[i_class]][-1,][,i_lambda]) # we remove the intercept
                    }
                }
            }
            for (i_lambda in 1:n_lambda) {
                for (i_class in 1:n_class) {
                    is.na(feature_coef_per_run[[i_lambda]][[i_class]]) <- feature_coef_per_run[[i_lambda]][[i_class]] == 0
                    feature_coef_per_lambda[[i_lambda]][[i_class]][,i_run] = rowMeans(feature_coef_per_run[[i_lambda]][[i_class]], na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                    feature_freq_per_lambda[[i_lambda]][[i_class]][,i_run] = rowSums(!is.na(feature_coef_per_run[[i_lambda]][[i_class]]))/n_fold
                }
            }
            pb$tick()
        }
        
        lambda_QF_est_all_runs = matrix(rep(NA,n_lambda*n_run),ncol=n_run)
        for (i_lambda in 1:n_lambda) {
            for (i_run in 1:n_run) {
                lambda_QF_est_all_runs[i_lambda,i_run] = QF(predicted_values_all_lambda[(1:n_instance)+(i_run-1)*n_instance,i_lambda],y)
            }
        }
        lambda_QF_est[[i_alpha]] = apply(lambda_QF_est_all_runs,1,median)
        best_lambda_index = which.max(lambda_QF_est[[i_alpha]])
        
        best_lambda[i_alpha] = lambda_values[[i_alpha]][best_lambda_index]
        model_QF_est[i_alpha] = lambda_QF_est[[i_alpha]][best_lambda_index]
        
        predicted_values[[i_alpha]] = matrix(rep(NA,n_instance*n_class),ncol=n_class)
        for (i_instance in 1:n_instance) {
            for (i_class in 1:n_class) {
                predicted_values[[i_alpha]][i_instance,i_class] =  sum(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
                best_lambda_index]==class[i_class])/n_run
            }
        }
        
        for (i_class in 1:n_class) {
            weight_tot = rowSums(feature_freq_per_lambda[[best_lambda_index]][[i_class]])
            feature_coef_wmean[[i_class]][,i_alpha] = rowSums(feature_coef_per_lambda[[best_lambda_index]][[i_class]]*feature_freq_per_lambda[[best_lambda_index]][[i_class]],na.rm=T)/weight_tot
            for (i_feature in 1:n_feature) {
                feature_coef_wsd[[i_class]][i_feature,i_alpha] = sqrt(sum(feature_freq_per_lambda[[best_lambda_index]][[i_class]][i_feature,]/weight_tot[i_feature] * (feature_coef_per_lambda[[best_lambda_index]][[i_class]][i_feature,] - feature_coef_wmean[[i_class]][i_feature,i_alpha])^2,na.rm=T))
            }
            feature_freq_mean[[i_class]][,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]][[i_class]],1,mean)
            feature_freq_sd[[i_class]][,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]][[i_class]],1,sd)
        }
        
        # permutation-based null
        null_QF_est_all_runs = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
        null_feature_coef = vector("list",n_class)
        null_feature_freq = vector("list",n_class)
        for (i_class in 1:n_class) {
            null_feature_coef[[i_class]] = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
            null_feature_freq[[i_class]] = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        }
        pb <- progress_bar$new(
        format = paste0("  running NULL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            for (i_perm_null in 1:n_perm_null) {   # permutation loop
                fold_distrib_OK = F
                fold_distrib_fail = 0
                while (!fold_distrib_OK) { # for categorical models, we must check that each class is in more than 1 fold.
                    y_RDM = sample(y) # we randomly permute the response vector
                    if (sum(apply(table(foldid_per_run[[i_alpha]][,i_run],y_RDM)>0,2,sum)>1)==n_class) {
                        fold_distrib_OK = T
                    } else {
                        fold_distrib_fail = fold_distrib_fail + 1
                        if (fold_distrib_fail>fold_distrib_fail.max) {
                            stop("Failure in class distribution across cross-validation folds. Increase n_fold, remove very small classes, or re-run with larger fold_distrib_fail.max\n")
                        }
                    }
                }
                null_predicted_values = rep(NA,n_instance)
                null_feature_coef_per_run = vector("list",n_class)
                for (i_fold in 1:n_fold) {
                    instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                    instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                    fit = glmnet(x[instance_in_bag,],y_RDM[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]][best_lambda_index])
                    null_predicted_values[which(instance_out_of_bag)] = predict(fit, x[instance_out_of_bag,], s=lambda_values[[i_alpha]][best_lambda_index], type=prediction_type)
                    for (i_class in 1:n_class) {
                        null_feature_coef_per_run[[i_class]] = cbind(null_feature_coef_per_run[[i_class]],coef(fit)[[i_class]][-1,]) # we remove the intercept
                    }
                }
                null_QF_est_all_runs[i_run,i_perm_null] = QF(null_predicted_values,y_RDM)
                for (i_class in 1:n_class) {
                    is.na(null_feature_coef_per_run[[i_class]]) <- null_feature_coef_per_run[[i_class]] == 0
                    null_feature_coef[[i_class]][,(i_run-1)*n_perm_null+i_perm_null] = rowMeans(null_feature_coef_per_run[[i_class]], na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                    null_feature_freq[[i_class]][,(i_run-1)*n_perm_null+i_perm_null] = rowSums(!is.na(null_feature_coef_per_run[[i_class]]))/n_fold
                }
            }
            pb$tick()
        }
        
        for (i_class in 1:n_class) {
            weight_tot = rowSums(null_feature_freq[[i_class]])
            null_feature_coef_wmean[[i_class]][,i_alpha] = rowSums(null_feature_coef[[i_class]]*null_feature_freq[[i_class]],na.rm=T)/weight_tot
            for (i_feature in 1:n_feature) {
                null_feature_coef_wsd[[i_class]][i_feature,i_alpha] = sqrt(sum(null_feature_freq[[i_class]][i_feature,]/weight_tot[i_feature] * (null_feature_coef[[i_class]][i_feature,] - null_feature_coef_wmean[[i_class]][i_feature,i_alpha])^2,na.rm=T))
            }
            null_feature_freq_mean[[i_class]][,i_alpha] = apply(null_feature_freq[[i_class]],1,mean)
            null_feature_freq_sd[[i_class]][,i_alpha] = apply(null_feature_freq[[i_class]],1,sd)
        }
        
        # Comparisons of model vs null
        n_tail = 0
        for (i_run in 1:n_run) {
            n_tail = n_tail + sum(null_QF_est_all_runs[i_run,]>lambda_QF_est_all_runs[best_lambda_index,i_run])
        }
        QF_model_vs_null_pval[i_alpha] = (n_tail+1)/(n_run*n_perm_null+1) # correction based on Phipson & Smyth (2010)
        
        for (i_class in 1:n_class) {
            for (i_feature in 1:n_feature) {
                n_coef = 0
                n_tail_coef = 0
                n_tail_freq = 0
                for (i_run in 1:n_run) {
                    coef_null_vs_model = abs(null_feature_coef[[i_class]][i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)])>=abs(feature_coef_per_lambda[[best_lambda_index]][[i_class]][i_feature,i_run])
                    n_coef = n_coef + sum(!is.na(coef_null_vs_model))
                    n_tail_coef = n_tail_coef + sum(coef_null_vs_model,na.rm=T)
                    n_tail_freq = n_tail_freq + sum(null_feature_freq[[i_class]][i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)]>=feature_freq_per_lambda[[best_lambda_index]][[i_class]][i_feature,i_run])
                }
                feature_coef_model_vs_null_pval[[i_class]][i_feature,i_alpha] = (n_tail_coef+1)/(n_coef+1)
                feature_freq_model_vs_null_pval[[i_class]][i_feature,i_alpha] = (n_tail_freq+1)/(n_run*n_perm_null+1)
            }
        }
        
        cat("alpha = ",alpha[i_alpha]," completed.\n")
    } # end of alpha loop
    
    # return object
    list(
    # input data and parameters
    predictor = x, response = y, alpha = alpha, family = family, nlambda = nlambda,
    nlambda.ext = nlambda.ext, seed = seed, scaled = scaled, n_fold = n_fold, n_run = n_run,
    n_perm_null = n_perm_null, QF_label = QF_label, instance = instance,
    feature = feature, fold_distrib_fail.max = fold_distrib_fail.max, glmnet_params = glmnet.control(),
    # summary results
    best_lambda = best_lambda, model_QF_est = model_QF_est, QF_model_vs_null_pval = QF_model_vs_null_pval,
    # detailed results for plots and downstream analysis
    lambda_values = lambda_values, lambda_QF_est = lambda_QF_est,
    predicted_values = predicted_values,
    feature_coef_wmean = as(feature_coef_wmean,"CsparseMatrix"), feature_coef_wsd = as(feature_coef_wsd,"CsparseMatrix"),
    feature_freq_mean = as(feature_freq_mean,"CsparseMatrix"), feature_freq_sd = as(feature_freq_sd,"CsparseMatrix"),
    null_feature_coef_wmean = as(null_feature_coef_wmean,"CsparseMatrix"), null_feature_coef_wsd = as(null_feature_coef_wsd,"CsparseMatrix"),
    null_feature_freq_mean = as(null_feature_freq_mean,"CsparseMatrix"), null_feature_freq_sd = as(null_feature_freq_sd,"CsparseMatrix"),
    feature_coef_model_vs_null_pval = as(feature_coef_model_vs_null_pval,"CsparseMatrix"),
    feature_freq_model_vs_null_pval = as(feature_freq_model_vs_null_pval,"CsparseMatrix")
    )
}

