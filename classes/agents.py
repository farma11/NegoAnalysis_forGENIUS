# coding: UTF-8
from xml.etree import ElementTree
import preference

class Agent(object):
    def __init__(self, agentID, name, prefPath):
        """Agentを生成する"""
        self.id = agentID
        self.name = name
        self.pref = preference.Preference(prefPath)
        self.acCNT = 0
        self.util = 0.0

        self.actions = []
        self.utilList = []

    def getID(self):
        """selfのAgentIDを返す"""
        return self.id

    def getName(self):
        """selfのnameを返す"""
        return self.name

    def getAcCNT(self):
        """selfのACの回数を返す"""
        return self.acCNT

    def incrementAcCNT(self):
        """selfのAC回数をインクリメントする"""
        self.acCNT += 1

    def setUtil(self, util):
        """selfの獲得効用値を入力"""
        self.util = util

    def getUtil(self):
        """selfの獲得効用値を返す"""
        return self.util

    def putActions(self, action):
        """selfの行動リストに行動を追加"""
        self.actions.append(action)

    def getActions(self):
        """selfの行動リストを返す"""
        return self.actions

    def getAction(self, rnd):
        """selfのrnd番目の行動を返す"""
        return self.actions[rnd-1]

    def calUtility(self, bid):
        """bidのUtilityを計算する"""
        return self.pref.calBidUtil(bid)

    def putUtilList(self, util):
        """selfの各Valueの効用値を格納"""
        self.utilList.append(util)

    def getUtilList(self):
        """selfの各Valueの効用値を返す"""
        return self.utilList

    def getPref(self):
        """selfの効用情報クラスを返す"""
        return self.pref


    def printResult(self):
        print "-- Agent " + str(self.id) + " --"
        print "AgentName  : " + str(self.name)
        print "AcceptN    : " + str(self.acCNT)
        print "Utility    : " + str(self.util)
