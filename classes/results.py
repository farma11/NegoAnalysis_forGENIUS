# coding: UTF-8
import agents, preference

class Result(object):
    def __init__(self, agentN):
        self.agentN = agentN # agent数
        self.agentList = []

        self.time = 0.0 # 合意形成にかかった時間
        self.rnd = 0 # 合意形成にかかったRound数
        self.pareto = -1.0 # 合意BidにおけるDist. pareto
        self.nash = -1.0 # 合意BidにおけるDist. nash
        self.sw = -1.0 # 合意形成Bidにおける社会的余剰

        self.bids = []

    def getAgentN(self):
        """selfにおけるAgent数を返す"""
        return self.agentN

    def setTime(self, time):
        """selfにおける交渉時間を格納"""
        self.time = time

    def getTime(self):
        """selfにおける交渉時間を返す"""
        return self.time

    def setRound(self, rnd):
        """selfにおける交渉Round数を格納"""
        self.rnd = rnd

    def getRound(self):
        """selfにおける交渉Round数を返す"""
        return self.rnd

    def setPareto(self, pareto):
        """selfにおける合意形成BidとのDist.Paretoを格納"""
        self.pareto = pareto

    def getPareto(self):
        """selfにおける合意形成BidとのDist.Paretoを返す"""
        return self.pareto

    def setNash(self, nash):
        """selfにおける合意形成BidとのDist.Nashを格納"""
        self.nash = nash

    def getNash(self):
        """selfにおける合意形成BidとのDist.Nashを格納"""
        return self.nash

    def setSocialWalfare(self, sw):
        """社会的余剰を格納"""
        self.sw = sw

    def getSocialWalfare(self):
        """社会的余剰を返す"""
        return self.sw

    def setAgentList(self, agentlist):
        """交渉参加者リストを格納"""
        self.agentList = agentlist

    def getAgentList(self):
        """交渉参加者リストを返す"""
        return self.agentList

    def getAgent(self, agentID):
        """agentIDを持つエージェントを返す"""
        return self.agentList[agentID-1]

    def putBids(self, bid):
        """selfの提案Bidを格納"""
        self.bids.append(bid)

    def getBids(self):
        """selfの提案Bid一覧（List）を返す"""
        return self.bids

    def getBid(self, idx):
        """"selfのインデックスidxの提案Bidを返す"""
        return self.bids[idx]

    def getBidSelectedAgent(self, agentID, rnd):
        """selfにおいて特定Agentの特定RoundのBidを返す"""
        return self.bids[(rnd - 1)*self.agentN + agentID - 1]

    def printResult(self):
        print "--- Negotiation Result ---"
        print "AgentN         : " + str(self.agentN)
        print "Time           : " + str(self.time)
        print "Round          : " + str(self.rnd)
        print "Dist. Pareto   : " + str(self.pareto)
        print "Dist. Nash     : " + str(self.nash)
        print "Social Walfare : " + str(self.sw)
        print "AcceptedBid    : " + str(self.bids[-1])
        for agent in self.agentList:
            agent.printResult()
            for opponent in self.agentList:
                print "Between " + opponent.getName() + " & me :[ABS] ",
                print preference.calDistPref_abs(agent.getPref(),opponent.getPref()),
                print "[SQ] ",
                print preference.calDistPref_sq(agent.getPref(),opponent.getPref())
