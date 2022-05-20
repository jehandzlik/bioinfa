########################################################
## R and R squared 
## 20-May-2022
########################################################

study_hours <- c(1,1,2,2,1,2,2,3,3,4,4,5)
current_grade <- c(65,78,76,76,79,80,81,84,88,85,96,90)
exam_score <- c(58,61,62,65,65,68,72,74,78,85,90,95)

df <- data.frame(StudyHours = study_hours, CurrentGrade = current_grade, ExamScore = exam_score)

########################################################
## Simple Linear Regression
plot(df$StudyHours, df$ExamScore)

# R or Pearson correlation coefficient
cor(df$StudyHours, df$ExamScore)
# R squared
cor(df$StudyHours, df$ExamScore)^2

lm_result <- lm(ExamScore ~ StudyHours, df)
summary(lm_result)

########################################################
## Multiple Linear Regression

## We can fit a multiple linear regression model using "study hours" and
## "current grade" as the predictor variables and "exam score" as the
## response variable.
multivariate_lm_results <- lm(ExamScore ~ StudyHours + CurrentGrade, df)
summary(multivariate_lm_results)

## R: Correlation between the actual exam scores and the predicted exam 
## scores made by the model is 0.978

cor(df$ExamScore, predict(multivariate_lm_results, df))

## R-squared: The R-squared for this regression model
## is 0.956. This tells us that 95.6% of the variation
## in the exam scores can be explained by the
## number of hours studied and the student's
## current grade in the class.
summary(multivariate_lm_results)

1-sum((df$ExamScore - predict(multivariate_lm_results, df))^2)/
  sum( (df$ExamScore - mean(df$ExamScore))^2 )
