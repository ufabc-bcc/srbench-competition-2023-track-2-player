import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn

Biks = pd.read_csv('datasets\\MontrealBikeLane.csv',index_col='Date',parse_dates=True)
weather = pd.read_csv('datasets\\WeatherInfo.csv',index_col='Date/Time',parse_dates=True)
print(Biks)
print(weather)
print(Biks.info())
print(weather.info())

Biks.isnull().sum()
# print(weather.info())


print(weather.describe())
Biks.plot(figsize=(15, 10))
plt.show()


#Differentiate between weekdays and weekends  0 is Monday
# berri_bikes = Biks
berri_bikes = Biks.iloc[:,:2].copy()
berri_bikes.index
berri_bikes.index.day
berri_bikes.index.weekday
berri_bikes.loc[:,'weekday'] = berri_bikes.index.weekday
print(berri_bikes)


all_df=berri_bikes.join(weather)
print(all_df)

all_df.groupby(['weekday'])['Berri1'].mean().plot(kind='line')
plt.show()

all_df.groupby(['Month'])['Berri1'].mean().plot(kind='line')
plt.show()

all_df.groupby(['Day'])['Berri1'].mean().plot(kind='line')
plt.show()

all_df.groupby(['Max Temp (°C)'])['Berri1'].mean().plot(kind='line')
plt.show()



corr = all_df.corr()
mask = np.array(corr)
mask[np.tril_indices_from(mask)] = False

plt.subplots(figsize=(10, 10))
sn.heatmap(corr, mask=mask,vmax=.8, square=True,annot=True)
plt.ylim(0,len(corr))
plt.tight_layout()

plt.show()
print(all_df.info())



from sklearn.model_selection import train_test_split

dropFeatures = ["Time","Year","Data Quality","Max Temp Flag","Min Temp Flag","Mean Temp (°C)","Mean Temp Flag","Heat Deg Days Flag",
                  "Cool Deg Days Flag","Total Rain (mm)","Total Rain Flag","Total Snow Flag","Total Precip (mm)","Total Precip Flag",
                 "Snow on Grnd (cm)","Snow on Grnd Flag","Dir of Max Gust (10s deg)","Dir of Max Gust Flag","Spd of Max Gust (km/h)","Spd of Max Gust Flag","Unnamed: 27"
               ,"Heat Deg Days (°C)","Cool Deg Days (°C)","Total Snow (cm)"
]

X1 = all_df.drop(dropFeatures, axis=1).copy()
X = X1.drop("Berri1", axis=1).copy()
yLabels = X1.drop(X, axis=1).copy()

train_x, test_x, train_y, test_y = train_test_split(X, yLabels, test_size=0.3, random_state=42)

train=pd.concat([train_x, train_y], axis=1)
test=pd.concat([test_x, test_y], axis=1)
train_x.to_csv('train.txt', index=False, sep='\t')
test_x.to_csv('test.txt', index=False, sep='\t')