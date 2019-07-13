# Project to analyze loan statistics and correlate profit margins
#	on loans to the characteristics of the individual who
#	requested the loan. 

# Read in the data 
loan_data = read.csv("LoanStats3a-organized.csv")
str(loan_data)

# When we take a look at the data, there are some catagories that
#	need some additional processing...
#	int_rate ==> convert to a float of %
rates = loan_data$int_rate
loan_data$int_rate = unlist(strsplit(as.vector(rates),"%"))
loan_data$int_rate = as.numeric(loan_data$int_rate)


# We now want to simplify the information in the data frame and only
#	keep the parameters and variables we plan to analyze
#	Variables to remove include...
#		funded_amnt ==> redundant
#		emp_title ==> too much variation in factors
#		issue_d ==> irrelevant 
#		desc ==> too much variation in factors
#		title ==> too much variation in factors
#		zip_code ==> incomplete information
#		earliest_cr_line ==> too many factors
#		total_pymnt_inv ==> redundant 
#		last_pymnt_d ==> irrelevant
#		revol_util ==> redundant 

# Create list of variables to remove
nonvars = c("funded_amnt", "emp_title", "issue_d", "desc", "title", "zip_code", "earliest_cr_line", "total_pymnt_inv", "last_pymnt_d", "revol_util")

# Reset the loan_data without those parameters
loan_data = loan_data[, !(names(loan_data) %in% nonvars) ]
str(loan_data)

# Here is a list of remaining categories and what they contain...
#
#	loan_amnt = size of the loan given (in $)
#	term = 'Factor' of "36 months" or "60 months" = length of the loan
#	int_rate = agreed upon interest rate for the loan (in %)
#	installment = agreed upon payment play (in $ per month)
#	grade = 'Factor' regarding the grading system of the lender
#	emp_length = 'Factor' of 12 levels indicating level of employment
#		= < 1 year, 1 year, ..., 9 years, 10+ years, n/a
#	home_ownership = 'Factor' of 5 levels indicating living accomodations
#	annual_inc = annual income of the individual
#	verification_status = 'Factor' of 3 levels indicating how the individual's
#		income was verified
#	loan_status = 'Factor' of 2 leves indicating whether or not the loan
#		was fully paid off
#	purpose = 'Factor' of 14 levels indicating the stated purpose for the loan
#	addr_state = 'Factor' of 50 levels indicating where the individual lives
#	dti = percentage of the individuals debts to income on a monthly basis
#	delinq_2yrs = number of times the individual has been delinq on payments
#	inq_last_6mths = number of credit inquires for the individual in last 6 month
#	open_acc = number of open credit lines for the individual
#	pub_rec = number of negative public records (crimes, bankruptcies, etc)
#	revol_bal = average credit card debt at end of billing cycle
#	total_acc = total number of credit lines the individual has had
#	total_pymnt = payment recieved from the individual for the loan
#	total_rec_prncp = payment recieved on the prinicple loan amount
#	total_rec_int = total amount of interest recieved on the loan
#	total_rec_late_fee = total amount of money recieved in late fees
#	last_pymnt_amnt = last monthly payment made by individual 


# For our analysis, we want only to consider the loans that have been paid off
#	To do this, we need to create a data subset on the $loan_status factor
paid_loans = subset(loan_data, loan_status == "Fully Paid")

# We can now remove loan_status from the dataframe
paid_loans = paid_loans[, !(names(paid_loans) == "loan_status") ]
str(paid_loans)

# Now we need to aggregate some data together for analysis
#	Create a new column called 'profit' and have this equal
#	the sum of the total_rec_int and total_rec_late_fee
paid_loans$profit = paid_loans$total_rec_int + paid_loans$total_rec_late_fee

# Then create a new column called 'roi' and compute
#	this as the ratio of the profit to the loan_amnt
paid_loans$roi = paid_loans$profit / paid_loans$loan_amnt


# Now we are ready for the analysis. Our goal here is to use the
#	roi and borrower characteristics to determine what type of
#	individual would yield the highest profit margins for lenders
#
# What characteristics should we consider? What columns correspond to 
#	information about the borrower?
#
#	emp_length, home_ownership, annual_inc, addr_state, dti, delinq_2yrs,
#	inq_last_6mths, open_acc, pub_rec, revol_bal, total_acc

# To test our analysis/model, we need to split out data set into a training
#	set of data and a test set of data.
# Use the 'sample.split' function and the 'dplyr' and 'caTools' library to split the
#	data such that 60% is used for training and 40% is used for testing
install.packages("dplyr")
library(dplyr)

install.packages("caTools")
library(caTools)

spl = sample.split(paid_loans$roi, SplitRatio = 0.6)
loanTrain = paid_loans %>% filter(spl == TRUE)
loanTest = paid_loans %>% filter(spl == FALSE)
str(loanTrain)

# Let's start by considering a linear regression model
LinMod = lm(roi ~ ., data=loanTrain)
summary(LinMod)

# Overall, we get an R^2 of 0.86 for this model, but many of the parameters
#	are insignificant. The parameters that are most correlated are...
#	
#	last_pymnt_amnt, total_acc, open_acc, dti, purpose, verification_status,
#	grade, loan_amnt, term, int_rate, installment
#
# Much of this should not be too surprising. If the int_rate is high and the
#	installments are low, then the roi will also be high. Likewise, if 
#	the lender put the individual themselves in the higher grade, then
#	they are confident that they will make profit. 

# What if we only consider the factors of importance that are also linked to
#	characteristics of the individual...
LinMod = lm(roi ~ total_acc + open_acc + dti + purpose, data=loanTrain)
summary(LinMod)

# The regression is no good. Perhaps a linear regression model is inappropriate
#	for this data set. 

# Next, we will make a regression tree using the 'rpart' libraries
install.packages("rpart")
library(rpart)
install.packages("rpart.plot")
library(rpart.plot)

loantree = rpart( roi ~ . , data = loanTrain)
prp(loantree)

# What this analysis reveals is that there really is almost no correlation
#	between the borrower's characteristics and the lenders roi
#	Instead, roi is based on the loan amount and the payment plan
#	set forth by the lender. 




# Let's consider now a different question, is there a relationship between the
#	borrower's characteristics and the interest rate given?
loantree = rpart( int_rate ~ . , data = loanTrain)
prp(loantree)

# A quick tree indicates that the lender's grade specification primarily controls
#	the interest rate. Confirm with a linear regression
LinMod = lm(int_rate ~ grade, data=loanTrain)
summary(LinMod)

# YES, the R^2 on the linear regression is 0.92 for only the grading factor
#	The higher the grade, the lower the interest rate 




# We can use this knowledge to our advantange now. How will an individual be
#	graded (i.e., what is their interest rate)
#	by the lender given their status and characteristics?
#
#	emp_length, home_ownership, annual_inc, addr_state, dti, delinq_2yrs,
#	inq_last_6mths, open_acc, pub_rec, revol_bal, total_acc

# By reductively removing factors and adding others, we are trying to determine the
#	parameter set most appropriate for predicting the interest rate 

LinMod = lm(int_rate ~ emp_length + loan_amnt + term + installment + purpose + verification_status + home_ownership + dti + delinq_2yrs + inq_last_6mths + total_acc + pub_rec + revol_bal, data=loanTrain)
summary(LinMod)

# The regression is not great (R^2 = 0.42), but this is the best result achievable using only
#	individual data and loan agreement information (thus, there is information on how someone is graded
#	that is not interpretable from the data)

loantree = rpart( int_rate ~ emp_length + loan_amnt + term + installment + purpose + verification_status + home_ownership + dti + delinq_2yrs + inq_last_6mths + total_acc + pub_rec + revol_bal, data = loanTrain)
prp(loantree)

# Based on the tree results, the most important factor in determining the grading
#	and/or interest rate to charge an individual are...
#
#	term, installment, total_acc, delinq_2yrs, and inq_last_6mths

# Summary: Individuals looking for short term loans with low installment plans, 
#	many prior credit lines, and few delinquinces and/or credit inquiries
#	will tend to be graded higher and get lower interest rates. 

# Now we can test the model using the predict function
predictTest = predict(loantree, newdata = loanTest)

# predictTest has now a number of observations for each of the bins in the tree
#	Each bin represents an approximate interest rate to charge
# Transform the numeric vector into a factor vector to view the numbers of 
#	unique factors in the prediction
factorizedResult = factor(predictTest)
summary(factorizedResult)

# Since the interest rate is correlated primarily with the grade, we will make a
#	confusion matrix of test grades vs. predicted interest rates
result = table(loanTest$grade, factorizedResult)
result

# This confusion matrix is not the best method of comparison though, but we do see
#	that the model does show that the majority of the low interest rate results
#	are predicted to be grade A. 

# Display the factor levels for our factorizedResult
fact_levels = levels(factorizedResult)
fact_levels


# Reform result as table of actual int_rate v factorizedResult
result = table(loanTest$int_rate, factorizedResult)
result

# Manually iterate through the table summing up all the correct results from
#	the ConfusionMatrix and storing that sum in x.
#	Assume that the int_rate (rows) of the table are to match the 
#	predictedFactorizedResult when the actual int_rates are less than
#	or equal to the prediction. 
j=0
i=0
x=0
for (val in result[,j+1]) { 
	cutoff = (as.numeric(colnames(result)[j+1]) + as.numeric(colnames(result)[j+2]))/2
	if (as.numeric(rownames(result)[i+1]) <= cutoff) { 
		#print(as.numeric(rownames(result)[i+1])) 
		x = x + val
	} 
	i = i+1 
}
x

j=1
i=0
for (val in result[,j+1]) { 
	cutoff1 = (as.numeric(colnames(result)[j]) + as.numeric(colnames(result)[j+1]))/2
	cutoff2 = (as.numeric(colnames(result)[j+1]) + as.numeric(colnames(result)[j+2]))/2
	if (as.numeric(rownames(result)[i+1]) > cutoff1) { 
	if (as.numeric(rownames(result)[i+1]) <= cutoff2) { 
		x = x + val
	} 
	}
	i = i+1 
}
x

j=2
i=0
for (val in result[,j+1]) { 
	cutoff1 = (as.numeric(colnames(result)[j]) + as.numeric(colnames(result)[j+1]))/2
	cutoff2 = (as.numeric(colnames(result)[j+1]) + as.numeric(colnames(result)[j+2]))/2
	if (as.numeric(rownames(result)[i+1]) > cutoff1) { 
	if (as.numeric(rownames(result)[i+1]) <= cutoff2) { 
		x = x + val
	} 
	}
	i = i+1 
}
x

j=3
i=0
for (val in result[,j+1]) { 
	cutoff1 = (as.numeric(colnames(result)[j]) + as.numeric(colnames(result)[j+1]))/2
	cutoff2 = (as.numeric(colnames(result)[j+1]) + as.numeric(colnames(result)[j+2]))/2
	if (as.numeric(rownames(result)[i+1]) > cutoff1) { 
	if (as.numeric(rownames(result)[i+1]) <= cutoff2) { 
		x = x + val
	} 
	}
	i = i+1 
}
x

j=4
i=0
for (val in result[,j+1]) { 
	cutoff1 = (as.numeric(colnames(result)[j]) + as.numeric(colnames(result)[j+1]))/2
	cutoff2 = (as.numeric(colnames(result)[j+1]) + as.numeric(colnames(result)[j+2]))/2
	if (as.numeric(rownames(result)[i+1]) > cutoff1) { 
	if (as.numeric(rownames(result)[i+1]) <= cutoff2) { 
		x = x + val
	} 
	}
	i = i+1 
}
x

j=5
i=0
for (val in result[,j+1]) { 
	cutoff = (as.numeric(colnames(result)[j]) + as.numeric(colnames(result)[j+1]))/2
	if (as.numeric(rownames(result)[i+1]) > cutoff) { 
		x = x + val
	} 
	i = i+1 
}
x

# Accuracy 
x/sum(result)

# Unfortunately, this model seems to be only about 38% accurate 
#	i.e., 62% of the time, the model does not predict the interest rate correctly 
