

import torch
from torch import nn
import torchvision.datasets as dsets
import torchvision.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.model_selection import train_test_split
import torch.utils.data as Data
import matplotlib.pyplot as plt
from torch.autograd import Variable
import torch.nn.functional as F
from torch.nn import init
import sys
sys.path.append('.../')
import resnet

#%%
x_dict = {}
path = '.../x'



for file in os.listdir(path):
    filename=os.fsdecode(file)
    x_dict[filename] = []
    f=open(f'{path}/{filename}')
    for line in f:
        s=line.strip().split('\t')
        x_dict[filename].append(list(np.array(s).astype(float)))
        
    f.close()

    
#%%
SS_dict={'C':[0,0,1],'H':[0,1,0], 'E':[1,0,0]}

x_ss = {}

for key, val in x_dict.items():
    ss = []
    for v in val:
        if v[:3] == [0.,0.,1.]:
            ss.append('C'+str(int(v[-1])))
        elif v[:3] == [0.,1.,0.]:
            ss.append('H'+str(int(v[-1])))
        else:
            ss.append('E'+str(int(v[-1])))
    x_ss[''.join(ss)] = key



#%%

y_dict = {}
path = '.../y'


for file in os.listdir(path):
    filename=os.fsdecode(file)
    y_dict[filename] = []
    f=open(f'{path}/{filename}')
    for line in f:
        s=line.strip().split('\t')
        y_dict[filename].append(list(np.array(s).astype(float)))
        
    f.close()
#%%

test_data = []
for key, val in x_dict.items():
    x = torch.tensor(val).float().view(1,-1)
    y = torch.tensor(y_dict[key]).float()
    test_data.append([x,y])
test_loader = torch.utils.data.DataLoader(dataset=test_data, batch_size=1, shuffle=True)
    

#%%
import math
from matplotlib.ticker import PercentFormatter
import matplotlib

confidence = {}


x_pred = {}



    

model = resnet.ResNet(in_channels=1)
model.load_state_dict(torch.load(f'.../parameters.pkl', map_location='cpu'),strict=True)

all_pred=[]
yt = []
bx_name = []
model.eval()
with torch.no_grad():
    for step, (b_x, b_y) in enumerate(test_loader):

        all_pred.append(model(b_x.type(torch.FloatTensor)).detach().numpy())
        yt.append(b_y)
        b_x = b_x.detach().numpy()
        ss = []
        for i in range(0,len(b_x[0][0]),4):
            s=list(b_x[0][0][i:i+4])
            if s[:3] == [0.,0.,1.]:
                ss.append('C'+str(int(s[-1])))
            elif s[:3] == [0.,1.,0.]:
                ss.append('H'+str(int(s[-1])))
            else:
                ss.append('E'+str(int(s[-1])))
        bx_name.append(x_ss[''.join(ss)])


def phi_correlation(a,b,c,d):
    return (a*b-c*d)/(((a+c)*(a+d)*(b+c)*(b+d))**0.5)

non_contact = [0]
average_contact_correlation = []
combine_name_correlation = {}
for i, e in enumerate(all_pred):
    a=0 #correct contact
    b=0 #correct non-contact
    c=0 #incorrectly assign to contact
    d=0 #incorrectly assign to non-contact
    for x,t in enumerate(e[0]):
        for m,n in enumerate(t):
            # asd[n]=0
            mn = np.float128(max(np.exp(n,dtype='float32')))
            if 1-mn <= 1e-7:
                conf = 7
            else:
                conf = int(-math.log10(1-mn))
            n=list(n)
            
            
            if yt[i][0][x][m] == n.index(max(n)):
                
                if conf in confidence:
                    confidence[conf].append(0)
                else:
                    confidence[conf] = [0]
                if n.index(max(n)) not in non_contact:
                    a+=1
                else:
                    b+=1
            else:
                if conf in confidence:
                    confidence[conf].append(1)
                else:
                    confidence[conf] = [1]
                if n.index(max(n)) not in non_contact:
                    c+=1
                else:
                    d+=1
    combine_name_correlation[bx_name[i]]=[phi_correlation(a,b,c,d)]
    average_contact_correlation.append(phi_correlation(a,b,c,d))

asd=np.mean(average_contact_correlation)

width_height_1 = (10, 5)
plt.figure(figsize=width_height_1)
plt.grid()
bins = np.arange(0, 1.1, 0.1)
matplotlib.rc('font', size=15)
plt.title(f'Average MCC {asd:.3f}')
matplotlib.rc('font', size=18)
plt.xlabel('Prediction accuracy (MCC)')
acc = sorted(average_contact_correlation)
plt.hist(average_contact_correlation,bins=bins, weights=np.ones(len(average_contact_correlation)) / len(average_contact_correlation))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

plt.xticks(np.arange(0,1.1,0.1))

plt.show()

contact_name_correlation = {}
non_contact = [0]
vas=2
correlation_x = {}
average_contact_correlation = []
contact = [1,2,3,4]
for i, e in enumerate(all_pred):
    a=0 #correct contact
    b=0 #correct non-contact
    c=0 #incorrectly assign to contact
    d=0 #incorrectly assign to non-contact
    for x,t in enumerate(e[0]):
        for m,n in enumerate(t):
            # asd[n]=0
            n=list(n)
            if yt[i][0][x][m] in contact and n.index(max(n)) in contact:
                

                a+=1
            elif yt[i][0][x][m] in non_contact and n.index(max(n)) in non_contact:
                b+=1
            elif yt[i][0][x][m] in contact and n.index(max(n)) in non_contact:
                d+=1
            else:
                c+=1
    contact_name_correlation[bx_name[i]] = [phi_correlation(a,b,c,d)]
    average_contact_correlation.append(phi_correlation(a,b,c,d))



asd=np.mean(average_contact_correlation)

width_height_1 = (10, 5)
plt.figure(figsize=width_height_1)
plt.grid()
bins = np.arange(0, 1.1, 0.1)
# bins.append(0.9)
matplotlib.rc('font', size=15)
plt.hist(average_contact_correlation,bins=bins, weights=np.ones(len(average_contact_correlation)) / len(average_contact_correlation))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

plt.title(f'Average MCC {asd:.3f}')
matplotlib.rc('font', size=18)
plt.xlabel('Prediction accuracy (MCC)')
acc = sorted(average_contact_correlation)
plt.xticks(np.arange(0,1.1,0.1))


plt.show()
    

    
    
#%%
combine_name_correlation={k: v for k, v in sorted(combine_name_correlation.items(), key=lambda item: item[1])}


contact_name_correlation={k: v for k, v in sorted(contact_name_correlation.items(), key=lambda item: item[1])}






#%%
percentage_ss = {}
f=open('.../percentage_ss')
for line in f:
    s=line.strip().split('\t')
    name = s[0]
    if name in combine_name_correlation:
        percentage_ss[name] = [float(s[1]),float(s[2])]

        combine_name_correlation[name].append(s[-1])
f.close()

#%%
max_sim = {}
f=open('.../secondary_similarity')
for line in f:
    s=line.strip().split('\t')
    
    max_sim[s[0]] = float(s[1])
f.close()
    
#%%
max_seq = {}
f=open('.../sequence_similarity')
for line in f:
    s=line.strip().split('\t')
    
    max_seq[s[0]] = float(s[1])
f.close()


#%%
x=[]
c1 = []
y1=[]
c=[]
y2 = []
y3 = []
for key, val in combine_name_correlation.items():
    # x.append(percentage_ss[key][0])
    c.append(percentage_ss[key][0])
    x.append(max_sim[key]) # secondary structure similarity
    # x.append(max_seq[key])
    y2.append(val[0])
    c1.append(percentage_ss[key][0])
    # y2.append(contact_name_correlation[key][0])
# x=c
x1 = []
x2 = []
for e,i in enumerate(y2):
    i=str(i)
    i1 = int(i[2])/10
    if i1 == 0:
        # print(i1)
        i1 = 1
    if float(i)<0.9:
        v=0
        for u in y2:
           if u>=i1 and u<i1+0.1:
               v+=1
        x1.append(v/len(y2))# percentage of MCC
    elif float(i)>=0.9 and float(i)<1:
        v=0
        for u in y2:
           if u>=i1 and u<=i1+0.1:
               v+=1
        x1.append(v/len(y2))# percentage of MCC
    elif float(i)==1:
        v=0
        for u in y2:
           if u>=i1-0.1 and u<=i1+0.1:
               v+=1
        x1.append(v/len(y2))# percentage of MCC

asd = {}
for i, e in enumerate(y2):
    qwe = int(e*10)
    if qwe == 10:
        qwe = 9
    if qwe in asd:
        asd[qwe].append(x[i])
    else:
        asd[qwe] = [x[i]]

for i, e in enumerate(y2):
    qwe = int(e*10)
    if qwe == 10:
        qwe = 9
    x2.append(np.mean(asd[qwe]))

width_height_1 = (10, 5)
plt.figure(figsize=width_height_1)
# plt.xticks([0.05*i for i in range(100)])
# plt.yticks([0.05*i for i in range(100)])

plt.grid()

# plt.ylabel('Percentage of correct predictions')
# plt.xlabel('confidence of predictions')
cmap = plt.cm.rainbow


# plt.title('Accuracy and confidence')
matplotlib.rc('font', size=18)
plt.ylabel('Fraction of accuracy')
matplotlib.rc('font', size=18)
plt.xlabel('Prediction accuracy (MCC)')
qwe = np.mean(y2)
# plt.plot(x1,y1,'.')
plt.scatter(y2,x1,c=np.array(x2),cmap=cmap, edgecolor='none',s=50*np.ones(len(y2)))
cbar = plt.colorbar()
cbar.set_label('Secondary structure similarity')
# cbar.set_label('Average maximum sequence identity')
# cbar.set_label('Average percentage of helices')

# cbar.set_label('Average percentage of helices')
plt.xticks(np.arange(0,1.1,0.1))
matplotlib.rc('font', size=15)
plt.title(f'Average MCC: {qwe:.3f}')

plt.show()

#%%
matplotlib.rc('font', size=18)
width_height_1 = (15, 8)
plt.figure(figsize=width_height_1)
plt.scatter(x,y2,c=c1,cmap=cmap, edgecolor='none',s=50*np.ones(len(y2)))

plt.xlabel('Maximum secondary structure similarity')
plt.ylabel('Prediction accuracy (MCC)')
cbar = plt.colorbar()
cbar.set_label('Average percentage of helices')
plt.grid()



#%%
xh = []
xe = []
xm = []
for key, val in combine_name_correlation.items():

    if percentage_ss[key][0]>percentage_ss[key][1]:
        xh.append(val[0])

    else:

        xe.append(val[0])

 
width_height_1 = (10, 5)
plt.figure(figsize=width_height_1)
plt.grid()
bins = np.arange(0, 1.1, 0.1)
# bins.append(0.9)

# plt.hist(xh,bins=bins, weights=np.ones(len(xh)) / len(xh))
plt.hist(xe,bins=bins, weights=np.ones(len(xe)) / len(xe))
# plt.hist(xm,bins=bins, weights=np.ones(len(xm)) / len(xm))

plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
matplotlib.rc('font', size=18)  
plt.xlabel('Prediction accuracy (MCC)')
asd = np.mean(xe)
# asd = np.mean(xm)
# asd = np.mean(xh)
matplotlib.rc('font', size=15)  
plt.title(f'Average MCC {asd:.3f}')
plt.xticks(np.arange(0,1.1,0.1))


plt.show()



#%%
matplotlib.rc('font', size=18)  

x = []
y = []
his = []
keys = sorted(list(confidence.keys()))
for key in keys:
    x.append(key)
    y.append(confidence[key].count(0)/len(confidence[key]))
    his+=list(key*np.ones(len(confidence[key])))
    
width_height_1 = (10, 5)
plt.figure(figsize=width_height_1)
plt.grid()
plt.plot(x,y,'-o')
plt.xlabel('Confidence')
plt.ylabel('Accuracy')

plt.show()
from matplotlib.ticker import PercentFormatter
width_height_1 = (10, 5)
plt.figure(figsize=width_height_1)
plt.grid()
plt.hist(his,weights=np.ones(len(his)) / len(his),bins=np.arange(0, 8, 0.999))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Confidence')
plt.ylabel('Distribution')





    
