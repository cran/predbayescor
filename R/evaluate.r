
## ratio.loss is the loss0to1 to loss1to0
evaluate_by_loss <- function( y.true, pred.prob, threshold=0.5)
{
     ratio.loss <- threshold / (1 - threshold)
     y.pred <- 1 * (pred.prob > threshold )

     loss1to0 <- (1+ratio.loss)/ratio.loss/2
     loss0to1 <- (1+ratio.loss)/2
     
     n1to0 <- sum(y.pred[y.true == 1] == 0)
     n0to1 <- sum(y.pred[y.true == 0] == 1)

     n <- length(y.true)

     loss <- (loss1to0 * n1to0 + loss0to1 * n0to1) / n

     sd <- sqrt((loss^2 * (n-n1to0-n0to1) +
                (loss1to0 - loss)^2 * n1to0 +
		(loss0to1 - loss)^2 * n0to1 )) / n

     attr(loss,"sd") <- sd
     
     loss
}
