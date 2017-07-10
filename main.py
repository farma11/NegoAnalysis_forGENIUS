# coding: UTF-8
import re

from xml.etree import ElementTree
from classes import agents, bids, results, preference, graphes

# 関数
def inputResultFloatData(line):
    target = re.compile("([0-9]+\.+[0-9]+)").search(line)
    if target != None:
        return float(target.group(1))
    else:
        return None

# 解析実行ファイルの読み取り
runfile = "./runfiles/test.xml"
runtree = ElementTree.parse(runfile)
runroot = runtree.getroot()

# Agent数, Log数の取得
agentN = int(runroot.find('agent_data').get('size'))
resultN = int(runroot.find('log').get('size'))

agentNameList = []
agentNames = runroot.find('agent_data')
for agent in agentNames:
    agentNameList.append(agent.get('name'))

# 効用情報Pathリストの初期化
prefPathList = []
prefPathRoot = runroot.find('preference').get('root')
for item in runroot.find('preference'):
    prefPathList.append(prefPathRoot + "/" + item.get('name'))

# ResultListの初期化
resultList = []
for i in range(1, resultN + 1):
    resultList.append(results.Result(agentN))

# 各解析用logについて
for logidx, log in enumerate(runroot.find('log')):
    # 解析用のlogのパスを取得, logファイルの読み込み
    logfile = runroot.find('log').get('root') + "/" + log.get('name')
    inputLog = open(logfile, 'r')

    # 現在の交渉LogについてAgentListを作成（データが揃ったのちにResultクラスに格納）
    nowAgentList = []
    for i in range(1, agentN + 1):
        nowAgentList.append(agents.Agent(i, agentNameList[i-1], prefPathList[i-1]))

    # Round毎の情報
    nowAgentID = 1
    nowRound = 1
    for line in inputLog:

        # Round xx
        roundN = re.compile("Round ([0-9]+)").search(line)
        if roundN != None:
            nowRound = int(roundN.group(1))
            resultList[logidx].setRound(nowRound)
            continue

        # Turn x: Party y ... START
        agentID = re.compile("Turn ([0-9]+):.*").search(line)
        if agentID != None:
            nowAgentID = int(agentID.group(1)) # 現在TurnのAgentIDの更新
            if line.find('(Accept)') != -1:
                nowAgentList[nowAgentID-1].putActions("Accept")
                nowAgentList[nowAgentID-1].incrementAcCNT()
                resultList[logidx].putBids(resultList[logidx].getBid(-1))
            elif line.find('(EndNegotiation)') != -1:
                nowAgentList[nowAgentID-1].putActions("EndNegotiation")
            else: # Offer
                nowAgentList[nowAgentID-1].putActions("Offer")
                bid = bids.divideValue(line)
                resultList[logidx].putBids(bid)
            # 現RoundのBid情報を格納
            for agent in nowAgentList:
                agent.putUtilList(agent.calUtility(resultList[logidx].getBid(-1)))

            continue

        # Result
        attribute = re.compile("(.+):").match(line)
        if attribute != None:
            att = attribute.group(1)
            if att == "Time (s)":
                resultList[logidx].setTime(inputResultFloatData(line))
            elif att == "Distance to pareto":
                resultList[logidx].setPareto(inputResultFloatData(line))
            elif att == "Distance to Nash":
                resultList[logidx].setNash(inputResultFloatData(line))
            elif att == "Social welfare":
                resultList[logidx].setSocialWalfare(inputResultFloatData(line))
            elif att == "Agent utility":
                util = inputResultFloatData(line)
                if util != None:
                    for agent in nowAgentList:
                        if agent.getUtil() == 0:
                            agent.setUtil(util)
                            break

    # 現在の交渉LogのAgentListをResultに格納
    resultList[logidx].setAgentList(nowAgentList)
    resultList[logidx].printResult()

    # グラフの描画
    graphes.plotNegoResult(resultList[logidx])

    inputLog.close()
