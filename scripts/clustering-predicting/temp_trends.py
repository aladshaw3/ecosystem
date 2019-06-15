## Python Data Analysis and Machine Learning Example ##
## Run python scripts using Python 3.5 or newer ##

''' Global Warming Analysis script:
    ----------------
    Object-Oriented approach to analyzing historical temperature data (1961 - 2009)
    and using that data to regress and train a model to predict future temperatures
    (2010 - 2016).
    
    Author:     Austin Ladshaw
    Date:       06/15/2019
    Copyright:  This software was designed and built by Austin Ladshaw.
    Copyright (c) 2019, all rights reserved.
'''

import pylab
import re
import random
import math

# cities in our weather data
CITIES = [
    'BOSTON',
    'SEATTLE',
    'SAN DIEGO',
    'PHILADELPHIA',
    'PHOENIX',
    'LAS VEGAS',
    'CHARLOTTE',
    'DALLAS',
    'BALTIMORE',
    'SAN JUAN',
    'LOS ANGELES',
    'MIAMI',
    'NEW ORLEANS',
    'ALBUQUERQUE',
    'PORTLAND',
    'SAN FRANCISCO',
    'TAMPA',
    'NEW YORK',
    'DETROIT',
    'ST LOUIS',
    'CHICAGO'
]

TRAINING_INTERVAL = range(1961, 2010)
TESTING_INTERVAL = range(2010, 2016)

class Climate(object):
    """
    The collection of temperature records loaded from given csv file
    """
    def __init__(self, filename):
        """
        Initialize a Climate instance, which stores the temperature records
        loaded from a given csv file specified by filename.

        Args:
            filename: name of the csv file (str)
        """
        self.rawdata = {}

        f = open(filename, 'r')
        header = f.readline().strip().split(',')
        for line in f:
            items = line.strip().split(',')

            date = re.match('(\d\d\d\d)(\d\d)(\d\d)', items[header.index('DATE')])
            year = int(date.group(1))
            month = int(date.group(2))
            day = int(date.group(3))

            city = items[header.index('CITY')]
            temperature = float(items[header.index('TEMP')])
            if city not in self.rawdata:
                self.rawdata[city] = {}
            if year not in self.rawdata[city]:
                self.rawdata[city][year] = {}
            if month not in self.rawdata[city][year]:
                self.rawdata[city][year][month] = {}
            self.rawdata[city][year][month][day] = temperature
            
        f.close()

    def random_day(self, leap=False):
        """
            Function to return two integers: random month and random day in that month
            Can include leap year if needed (Defaults to False)
            
            Args:
                leap: False = neglect leap years (Excludes possibility of returning Feb. 29)
                
            Returns:
                a tuple containing an int for the month and an int for the day of the month
        """
        month = random.randint(1,12)
        day = 0
        if month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12:
            day = random.randint(1,31)
        elif month == 2:
            if leap == False:
                day = random.randint(1,28)
            else:
                day = random.randint(1,29)
        else:
            day = random.randint(1,30)
        return (month, day)

    def get_yearly_temp(self, city, year):
        """
        Get the daily temperatures for the given year and city.

        Args:
            city: city name (str)
            year: the year to get the data for (int)

        Returns:
            a 1-d pylab array of daily temperatures for the specified year and
            city
        """
        temperatures = []
        assert city in self.rawdata, "provided city is not available"
        assert year in self.rawdata[city], "provided year is not available"
        for month in range(1, 13):
            for day in range(1, 32):
                if day in self.rawdata[city][year][month]:
                    temperatures.append(self.rawdata[city][year][month][day])
        return pylab.array(temperatures)

    def get_daily_temp(self, city, month, day, year):
        """
        Get the daily temperature for the given city and time (year + date).

        Args:
            city: city name (str)
            month: the month to get the data for (int, where January = 1,
                December = 12)
            day: the day to get the data for (int, where 1st day of month = 1)
            year: the year to get the data for (int)

        Returns:
            a float of the daily temperature for the specified time (year +
            date) and city
        """
        assert city in self.rawdata, "provided city is not available"
        assert year in self.rawdata[city], "provided year is not available"
        assert month in self.rawdata[city][year], "provided month is not available"
        assert day in self.rawdata[city][year][month], "provided day is not available"
        return self.rawdata[city][year][month][day]

def se_over_slope(x, y, estimated, model):
    """
    For a linear regression model, calculate the ratio of the standard error of
    this fitted curve's slope to the slope. The larger the absolute value of
    this ratio is, the more likely we have the upward/downward trend in this
    fitted curve by chance.
    
    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        estimated: an 1-d pylab array of values estimated by a linear
            regression model
        model: a pylab array storing the coefficients of a linear regression
            model

    Returns:
        a float for the ratio of standard error of slope to slope
    """
    assert len(y) == len(estimated)
    assert len(x) == len(estimated)
    EE = ((estimated - y)**2).sum()
    var_x = ((x - x.mean())**2).sum()
    SE = pylab.sqrt(EE/(len(x)-2)/var_x)
    return SE/model[0]

def isclose(a,b):
    if abs(a - b) < 1e-6:
        return True
    else:
        return False

def generate_models(x, y, degs):
    """
    Generate regression models by fitting a polynomial for each degree in degs
    to points (x, y).

    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        degs: a list of degrees of the fitting polynomial

    Returns:
        a list of pylab arrays, where each array is a 1-d array of coefficients
        that minimizes the squared error of the fitting polynomial
    """
    models = []
    try:
        for d in degs:
            models.append(pylab.polyfit(x,y,d))
    except:
        models.append(pylab.polyfit(x,y,degs))
    return models

def r_squared(y, estimated):
    """
    Calculate the R-squared error term.
    
    Args:
        y: 1-d pylab array with length N, representing the y-coordinates of the
            N sample points
        estimated: an 1-d pylab array of values estimated by the regression
            model

    Returns:
        a float for the R-squared error term
    """
    mean = pylab.mean(y)
    return 1.0 - (pylab.sum((y-estimated)**2) / pylab.sum((y - mean)**2))

def evaluate_models_on_training(x, y, models):
    """
    For each regression model, compute the R-squared value for this model with the
    standard error over slope of a linear regression line (only if the model is
    linear), and plot the data along with the best fit curve.

    For the plots, you should plot data points (x,y) as blue dots and your best
    fit curve (aka model) as a red solid line. You should also label the axes
    of this figure appropriately and have a title reporting the following
    information:
        degree of your regression model,
        R-square of your model evaluated on the given data points,
        and SE/slope (if degree of this model is 1 -- see se_over_slope). 

    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        models: a list containing the regression models you want to apply to
            your data. Each model is a pylab array storing the coefficients of
            a polynomial.

    Returns:
        None
    """
    '''
    poly = pylab.poly1d(models[0])
    print poly
    f_vals = []
    for i in range(0,len(x)):
        f_vals.append(poly(x[i]))
    return f_vals
    '''
    poly = []
    f_vals = []
    m = 0
    for model in models:
        poly.append(pylab.poly1d(model))
        f_vals.append([])
        for i in range(0,len(x)):
            f_vals[m].append(poly[m](x[i]))
        pylab.figure()
        pylab.plot(x,y,'bo',label = 'Data')
        pylab.plot(x,f_vals[m],'r', label = 'Model')
        pylab.legend(loc='upper left')
        pylab.xlabel("Years")
        pylab.ylabel("Degrees Celsius")
        axes = pylab.gca()
        axes.set_xlim([x.min()-1.0,x.max()+1.0])
        axes.set_ylim([y.min()-1.0,y.max()+1.0])
        deg = len(model)-1
        r2 = r_squared(y,f_vals[m])
        if deg > 1:
            pylab.title("Polynomial Degree = " + str(deg) + "\n" + "R^2 = " + str(r2))
        else:
            se = se_over_slope(x, y, f_vals[m], model)
            sig = False
            if abs(se) < 0.5:
                sig = True
            pylab.title("Polynomial Degree = " + str(deg) + "\n" + "R^2 = " + str(r2) + "\nSE = " + str(se))
        m += 1
    pylab.show()

def gen_cities_avg(climate, multi_cities, years):
    """
    Compute the average annual temperature over multiple cities.

    Args:
        climate: instance of Climate
        multi_cities: the names of cities we want to average over (list of str)
        years: the range of years of the yearly averaged temperature (list of
            int)

    Returns:
        a pylab 1-d array of floats with length = len(years). Each element in
        this array corresponds to the average annual temperature over the given
        cities for a given year.
    """
    national_avg = []
    for year in years:
        city_avg = 0.0
        for city in multi_cities:
            city_avg += pylab.mean(climate.get_yearly_temp(city,year))
        national_avg.append(city_avg/float(len(multi_cities)))
    return national_avg

def moving_average(y, window_length):
    """
    Compute the moving average of y with specified window length.

    Args:
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        window_length: an integer indicating the window length for computing
            moving average

    Returns:
        an 1-d pylab array with the same length as y storing moving average of
        y-coordinates of the N sample points
    """
    avg = []
    for i in range(0,len(y)):
        #first window_length - 1 elements are scoped inwards
        if i < (window_length-1):
            avg.append( pylab.mean(y[0:i+1]) )
        else:
            avg.append( pylab.mean(y[i-window_length+1:i+1]) )
    return avg

def rmse(y, estimated):
    """
    Calculate the root mean square error term.

    Args:
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        estimated: an 1-d pylab array of values estimated by the regression
            model

    Returns:
        a float for the root mean square error term
    """
    return math.sqrt(pylab.sum((y-estimated)**2) / float(len(y)))

def gen_std_devs(climate, multi_cities, years):
    """
    For each year in years, compute the standard deviation over the averaged yearly
    temperatures for each city in multi_cities. 

    Args:
        climate: instance of Climate
        multi_cities: the names of cities we want to use in our std dev calculation (list of str)
        years: the range of years to calculate standard deviation for (list of int)

    Returns:
        a pylab 1-d array of floats with length = len(years). Each element in
        this array corresponds to the standard deviation of the average annual 
        city temperatures for the given cities in a given year.
    """
    '''
    national_std = []
    for year in years:
        city_avg = 0.0
        for city in multi_cities:
            city_avg += pylab.mean(climate.get_yearly_temp(city,year))
        national_avg = city_avg/float(len(multi_cities))
        
        sqsum = 0.0
        count = 0
        for city in multi_cities:
            for temp in climate.get_yearly_temp(city,year):
                sqsum += (temp - national_avg)**2
                count += 1
        
        national_std.append(math.sqrt(sqsum/float(count)))
        
    return national_std
    '''
    national_std = []
    y = 0
    daily_list = []
    for year in years:
        daily_list.append([])
        for month in range(1,13):
            for day in range(1, 32):
                daily = []
                for city in multi_cities:
                    try:
                        daily.append(climate.get_daily_temp(city, month, day, year))
                    except:
                        pass
                if len(daily) != 0:
                    avg_daily = pylab.mean(pylab.array(daily))
                    daily_list[y].append(avg_daily)
        national_std.append(pylab.std(pylab.array(daily_list[y])))
        y += 1
    return national_std

def evaluate_models_on_testing(x, y, models):
    """
    For each regression model, compute the RMSE for this model and plot the
    test data along with the model’s estimation.

    For the plots, you should plot data points (x,y) as blue dots and your best
    fit curve (aka model) as a red solid line. You should also label the axes
    of this figure appropriately and have a title reporting the following
    information:
        degree of your regression model,
        RMSE of your model evaluated on the given data points. 

    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        models: a list containing the regression models you want to apply to
            your data. Each model is a pylab array storing the coefficients of
            a polynomial.

    Returns:
        None
    """
    poly = []
    f_vals = []
    m = 0
    for model in models:
        poly.append(pylab.poly1d(model))
        f_vals.append([])
        for i in range(0,len(x)):
            f_vals[m].append(poly[m](x[i]))
        pylab.figure()
        pylab.plot(x,y,'bo',label = 'Data')
        pylab.plot(x,f_vals[m],'r', label = 'Model')
        pylab.legend(loc='upper left')
        pylab.xlabel("Years")
        pylab.ylabel("Degrees Celsius")
        axes = pylab.gca()
        axes.set_xlim([x.min()-1.0,x.max()+1.0])
        axes.set_ylim([y.min()-1.0,y.max()+1.0])
        deg = len(model)-1
        rm = rmse(y,f_vals[m])
        pylab.title("Prediction:\nPolynomial Degree = " + str(deg) + "\n" + "RSME = " + str(rm))
        m += 1
    pylab.show()

if __name__ == '__main__':

    # Part A.4 - Data Training Example
    '''
    print("Part A.4-1 : Temperature of New York on Jan. 10th from 1961 to 2009")
    climate_data = Climate("data.csv")
    month = 1
    day = 10
    dates = TRAINING_INTERVAL
    temps = []
    for year in dates:
        temps.append(climate_data.get_daily_temp('NEW YORK',month,day,year))
    model = generate_models(pylab.array(dates), pylab.array(temps), 1)
    evaluate_models_on_training(pylab.array(dates), pylab.array(temps),model)
    print()

    print("Part A.4-2 : Average Temperature of New York from 1961 to 2009")
    avg_temps = []
    for year in dates:
        avg_temps.append( pylab.mean(climate_data.get_yearly_temp('NEW YORK',year)) )
    model2 = generate_models(pylab.array(dates), pylab.array(avg_temps), 1)
    evaluate_models_on_training(pylab.array(dates), pylab.array(avg_temps),model2)
    print()
    '''

    # Part B
    '''
    print("Part B : National Average Temperature from 1961 to 2009")
    climate_data = Climate("data.csv")
    city_list = CITIES
    dates = TRAINING_INTERVAL
    national_avg = gen_cities_avg(climate_data,city_list,dates)
    model3 = generate_models(pylab.array(dates), pylab.array(national_avg), 1)
    evaluate_models_on_training(pylab.array(dates), pylab.array(national_avg),model3)
    print()
    '''
    
    # Part C
    '''
    print("Part C : National 5-year Average Temperature from 1961 to 2009")
    climate_data = Climate("data.csv")
    city_list = CITIES
    dates = TRAINING_INTERVAL
    national_avg = gen_cities_avg(climate_data,city_list,dates)
    five_year_avg = moving_average(national_avg,5)
    model4 = generate_models(pylab.array(dates), pylab.array(five_year_avg), 1)
    evaluate_models_on_training(pylab.array(dates), pylab.array(five_year_avg),model4)
    print(list)
    print()
    '''

    # Part D.2
    print("Part D : Predict Temperatures for 2010 to 2015 using Training Data")
    climate_data = Climate("data.csv")
    city_list = CITIES
    training_years = pylab.array(TRAINING_INTERVAL)
    testing_years = pylab.array(TESTING_INTERVAL)
    
    #Training Models
    national_avg = gen_cities_avg(climate_data,city_list,training_years)
    five_year_avg = moving_average(national_avg,5)
    model5 = generate_models(pylab.array(training_years), pylab.array(five_year_avg), [1,2,20])
    evaluate_models_on_training(pylab.array(training_years), pylab.array(five_year_avg),model5)
    
    #Testing Models
    future_avg = gen_cities_avg(climate_data,city_list,testing_years)
    future_5_year = moving_average(future_avg,5)
    evaluate_models_on_testing(pylab.array(testing_years), pylab.array(future_5_year),model5)
    print()

    # Part E
    print("Part E : Show Trend in Increasing Extreme Weather from 1961 to 2009")
    national_std = gen_std_devs(climate_data,city_list,training_years)
    std_5_yr = moving_average(national_std,5)
    extreme = pylab.array(std_5_yr) + pylab.array(five_year_avg)
    model6 = generate_models(pylab.array(training_years), pylab.array(extreme), 1)
    evaluate_models_on_training(pylab.array(training_years), pylab.array(extreme),model6)
    print()

