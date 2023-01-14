#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 17:38:45 2022

@author: runfeng
"""

import torch
from torch import nn
import torchvision.datasets as dsets
import torchvision.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.model_selection import train_test_split
import torch.utils.data as Data
from torch.autograd import Variable
import torch.nn.functional as F
from torch.nn import init


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


class Conv3x3(torch.nn.Module):
    def __init__(self,In_channel,Med_channel,Out_channel,downsample=False):
        super(Conv3x3, self).__init__()
        self.stride = 1
        if downsample == True:
            self.stride = 2
            # self.stride = 3


        self.layer = torch.nn.Sequential(
            torch.nn.Conv1d(In_channel, Med_channel, 3,stride=self.stride,padding=1),
            torch.nn.BatchNorm1d(Med_channel),
            torch.nn.ReLU(),
            torch.nn.Conv1d(Med_channel, Med_channel, 3,stride=1,padding=1),
            torch.nn.BatchNorm1d(Med_channel),
            torch.nn.ReLU(),
            torch.nn.Conv1d(Med_channel, Out_channel, 3,stride=1,padding=1),
            torch.nn.BatchNorm1d(Out_channel),
            torch.nn.ReLU(),
        )

        if In_channel != Out_channel:
            self.res_layer = torch.nn.Conv1d(In_channel, Out_channel,3,stride=self.stride,padding=1)
        else:
            self.res_layer = None

    def forward(self,x):
        # print(x.shape)
        if self.res_layer is not None:
            residual = self.res_layer(x)
        else:
            residual = x
        return self.layer(x)+residual
    
    
    
#%%

class ResNet(torch.nn.Module):
    def __init__(self,in_channels):
        super(ResNet, self).__init__()

        self.features = torch.nn.Sequential(
            torch.nn.Conv1d(in_channels,64,kernel_size=3,stride=1,padding=1),
            torch.nn.BatchNorm1d(64),
            torch.nn.ReLU(),
            torch.nn.MaxPool1d(kernel_size=3,stride=1,padding=1),
            
            
            
            Conv3x3(64,64,256,False),
            Conv3x3(256,64,256,False),
            Conv3x3(256,64,256,False),
            #
            Conv3x3(256,128,512, True),
            Conv3x3(512,128,512, False),
            Conv3x3(512,128,512, False),
            #
            Conv3x3(512,256,1024, True),
            Conv3x3(1024,256,1024, False),
            Conv3x3(1024,256,1024, False),
            #
            Conv3x3(1024,512,2048, True),
            Conv3x3(2048,512,2048, False),
            Conv3x3(2048,512,2048, False),
            
            #
            Conv3x3(2048,1024,4096, False),
            Conv3x3(4096,1024,4096, False),
            Conv3x3(4096,1024,4096, False),


            torch.nn.AvgPool1d(kernel_size=3,stride=1,padding=1)
        )
        self.classifer = torch.nn.Sequential(

            torch.nn.Conv2d(4096,5,kernel_size=3,stride=1,padding=1),

        )
       

    def forward(self,x):

        x = self.features(x)

        x = x[:,:,:,None]
        x = torch.matmul(x,torch.transpose(x,-1,-2))


        x = self.classifer(x)

        x = torch.transpose(x,-1,-3)
        x = torch.log_softmax(x,dim=-1)

        return x



