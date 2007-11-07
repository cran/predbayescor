prpr <- function(x,d)
{  if(!is.nan(x)) 
      sprintf(paste("%1.",d,"f",sep=""),x)
   else sprintf(paste("%",d+2,"s",sep=""),"--")
}
present.pair.tex <- function(tc1,tu1,no.sel1,tc2,tu2,no.sel2,no.all,
                             ttc1,ttu1,ttc2,ttu2)
{  
    cat("\\begin{tabular}\n");
    cat("{c|ccc|ccc||ccc|ccc}\n")
    
    cat("\\hline\n")
    if(no.sel1==1) fth1 = "feature" else fth1="features"
    if(no.sel2==1) fth2 = "feature" else fth2="features"

    cat("        ","&",
        "\\multicolumn{6}{c}{",no.sel1,fth1,"selected out of",no.all,"}",
	           "&",
	"\\multicolumn{6}{c}{",no.sel2,fth2,"selected out of",no.all,"}",
	"\\\\ \n")
    cat("\\hline\n")
    cat("        ","&","\\multicolumn{3}{c}{Corrected}",
                   "&","\\multicolumn{3}{c}{Uncorrected}",
	           "&","\\multicolumn{3}{c}{Corrected}",
                   "&","\\multicolumn{3}{c}{Uncorrected}",
	"\\\\ \n");

    cat("\\hline\n")
    cat("TIME",    "&","\\multicolumn{3}{c}{",prpr(ttc1,0),"}",
                   "&","\\multicolumn{3}{c}{",prpr(ttu1,0),"}",
		   "&","\\multicolumn{3}{c}{",prpr(ttc2,0),"}",
		   "&","\\multicolumn{3}{c}{",prpr(ttu2,0),"}",
        "\\\\ \n");
    cat("\\hline\n")

    cat("AMLL",    "&","\\multicolumn{3}{c}{",prpr(tc1$aml,3),"}",
                   "&","\\multicolumn{3}{c}{",prpr(tu1$aml,3),"}",
		   "&","\\multicolumn{3}{c}{",prpr(tc2$aml,3),"}",
		   "&","\\multicolumn{3}{c}{",prpr(tu2$aml,3),"}",
        "\\\\ \n");
    cat("\\hline\n")

    cat("MSE","&","\\multicolumn{3}{c}{",prpr(tc1$mse,3),"}",
                   "&","\\multicolumn{3}{c}{",prpr(tu1$mse,3),"}",
		   "&","\\multicolumn{3}{c}{",prpr(tc2$mse,3),"}",
		   "&","\\multicolumn{3}{c}{",prpr(tu2$mse,3),"}",
        "\\\\ \n");
    cat("\\hline\n")
    
    cat("CATEGORY","&","NO.","&","\\small{Pred.}","&","\\small{Frac.1}",
                   "&","NO.","&","\\small{Pred.}","&","\\small{Frac.1}",
	           "&","NO.","&","\\small{Pred.}","&","\\small{Frac.1}",
	           "&","NO.","&","\\small{Pred.}","&","\\small{Frac.1}",
	"\\\\ \n");
    
    for(i in 1:10)
    cat(sprintf("%1.1f",(i-1)/10),"-",sprintf("%1.1f",i/10),  
        "&",sprintf("%4.0f",tc1$summary.pred[i,2]), "&", prpr(tc1$summary.pred[i,3],3), "&",prpr(tc1$summary.pred[i,4],3),
        "&",sprintf("%4.0f",tu1$summary.pred[i,2]), "&", prpr(tu1$summary.pred[i,3],3), "&",prpr(tu1$summary.pred[i,4],3),
        "&",sprintf("%4.0f",tc2$summary.pred[i,2]), "&", prpr(tc2$summary.pred[i,3],3), "&",prpr(tc2$summary.pred[i,4],3),
        "&",sprintf("%4.0f",tu2$summary.pred[i,2]), "&", prpr(tu2$summary.pred[i,3],3), "&",prpr(tu2$summary.pred[i,4],3),
"\\\\ \n");
    cat("\\hline\n")

    cat("\\end{tabular}\n");
   
    
}

present.sgl.tex <- function(tc1,no.all,ttc1)
{   
    cat("\\begin{tabular}[t]\n");
    cat("{c|ccc}\n")
    
    cat("\\hline\n")
    
    cat("        ","&",
        "\\multicolumn{3}{c}{Complete data}",
	"\\\\ \n")
    cat("\\hline\n")
    cat("\\\\ \n");
    cat("\\hline\n")
    cat("TIME",    "&","\\multicolumn{3}{c}{",prpr(ttc1,0),"}",
        "\\\\ \n");

    cat("\\hline\n")    
    cat("AMLL",    "&","\\multicolumn{3}{c}{",prpr(tc1$aml,3),"}",
        "\\\\ \n");
    cat("\\hline\n")
    cat("MSE","&","\\multicolumn{3}{c}{",prpr(tc1$mse,3),"}",
        "\\\\ \n");
    cat("\\hline\n")

    cat("CATEGORY","&","NO.","&","\\small{Pred.}","&","\\small{Frac.1}",
                   	"\\\\ \n");
    

    for(i in 1:10)
    cat(sprintf("%1.1f",(i-1)/10),"-",sprintf("%1.1f",i/10),  
        "&",sprintf("%4.0f",tc1$summary.pred[i,2]), "&", prpr(tc1$summary.pred[i,3],3), "&",prpr(tc1$summary.pred[i,4],3),
                   "\\\\ \n");
    cat("\\hline\n")

    cat("\\end{tabular}\n");
   
    
}
