# coding: UTF-8
import results, agents, preference, bids

# グラフ化に必要なものの準備
import pylab
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy import genfromtxt
from matplotlib import cm

# データの扱いに必要なライブラリ
import pandas as pd
import numpy as np
import datetime as dt

import preference

# pandas用
plt.style.use('ggplot')
font = {'family' : 'meiryo'}
matplotlib.rc('font', **font)

# 交渉結果を描画する
def plotNegoResult(res):
    agentList = res.getAgentList()

    util = [
        agentList[0].getUtilList(),
        agentList[1].getUtilList(),
        agentList[2].getUtilList()
    ]
    colors = ["r","g","b"]

    gs = gridspec.GridSpec(2,3)
    ax = [
        plt.subplot(gs[0,0]), plt.subplot(gs[0,1]),plt.subplot(gs[0,2]), # 2者間
        plt.subplot(gs[1,0:2]),                                          # 3者間折れ戦
        plt.subplot(gs[1,2], projection="3d")                            # 3D散布図
    ]

    # ２者間での推移のグラフ化
    for i in range(3):
        j = (i+1)%3
        ax[i].plot(util[i][i::3], util[j][i::3], marker=".", color=colors[i], label=str(agentList[i].getName()))
        ax[i].plot(util[i][j::3], util[j][j::3], marker=".", color=colors[j], label=str(agentList[j].getName()))
        #ax[i].set_title("PrefDist: " + str(preference.calDistPref_sq(agentList[i].getPref(),agentList[j].getPref())))
        ax[i].set_title(str(res.getAgent(i+1).getName()) + " vs " + str(res.getAgent(j+1).getName()))
        ax[i].set_xlabel(str(agentList[i].getName()) + "\'s Utility")
        ax[i].set_ylabel(str(agentList[j].getName()) + "\'s Utility")
        ax[i].set_ylim(-0.05, 1.05)
        ax[i].set_xlim(-0.05, 1.05)
        ax[i].legend(loc="lower left")
        ax[i].grid(True)

    # 3者間での効用値の推移をグラフ化
    totalT = res.getTime()
    t = [totalT * i / (len(util[0])+1) for i in range(1, len(util[0])+1)]
    ax[3].plot(t, util[0], marker=".", color=colors[0], label=str(agentList[0].getName()))
    ax[3].plot(t, util[1], marker=".", color=colors[1], label=str(agentList[1].getName()))
    ax[3].plot(t, util[2], marker=".", color=colors[2], label=str(agentList[2].getName()))
    ax[3].set_title("3-party Negotiation")
    ax[3].set_xlabel('time')
    ax[3].set_ylabel('utility')
    ax[3].set_ylim(-0.05, 1.05)
    ax[3].set_xlim(0.0, totalT)
    ax[3].legend(loc="lower right")
    ax[3].grid(True)

    # 3D散布図の描画
    ax[4].scatter3D(util[0][0::3], util[1][0::3], util[2][0::3], marker="x", color=colors[0], label=str(agentList[0].getName()))
    ax[4].scatter3D(util[0][1::3], util[1][1::3], util[2][1::3], marker="o", color=colors[1], label=str(agentList[1].getName()))
    ax[4].scatter3D(util[0][2::3], util[1][1::3], util[2][1::3], marker="^", color=colors[2], label=str(agentList[2].getName()))
    ax[4].set_xlabel(str(agentList[0].getName()))
    ax[4].set_ylabel(str(agentList[1].getName()))
    ax[4].set_zlabel(str(agentList[2].getName()))
    ax[4].set_xlim(-0.05, 1.05)
    ax[4].set_ylim(-0.05, 1.05)
    ax[4].set_zlim(-0.05, 1.05)
    ax[4].legend(loc="upper left")
    ax[4].grid(True)

    plt.show()
