evaluate_by_loss <- function( y.true, pred.prob, ratio.loss=10)
{

     threshold <- 1/(1+ratio.loss)

     y.pred <- 1 * (pred.prob > threshold )

     if( ratio.loss <= 1) {
         loss1to0 <- 1
         loss0to1 <- 1/ratio.loss
     }
     else
     {
	 loss1to0 <- ratio.loss
         loss0to1 <- 1

     }
     n1to0 <- sum(y.pred[y.true == 1] == 0)

     n0to1 <- sum(y.pred[y.true == 0] == 1)

     n <- length(y.true)

     loss <- (loss1to0 * n1to0 + loss0to1 * n0to1)/ n

     sd <- sqrt((loss^2 * (n-n1to0-n0to1) +
                (loss1to0 - loss)^2 * n1to0 +
		(loss0to1 - loss)^2 * n0to1 )) / n

     list(sd=sd,loss=loss)


}
