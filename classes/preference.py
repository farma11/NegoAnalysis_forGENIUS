# coding: UTF-8
from xml.etree import ElementTree
import multiprocessing as mp # 並列化
import sys
import copy
import decimal

class Preference(object):
    def __init__(self, xmlfile):
        self.xml = xmlfile

        self.tree = ElementTree.parse(self.xml)
        self.treeRoot = self.tree.getroot()
        self.utils = self.calValueUtils()
        self.allBidSize = self.calAllBidSize()

    def getTreeRoot(self):
        """xmlの木のrootを返す"""
        return self.treeRoot

    def getIssueSize(self):
        """Issueの数を返す"""
        return len(self.treeRoot.find('objective').findall('issue'))

    def getValueSize(self, issueID):
        """あるIssueのValueの数を返す"""
        return len(self.treeRoot.find('objective').findall('issue')[issueID-1].findall('item'))


    def getIssueWeight(self):
        """各Issueの重みリストを返す"""
        weights = []
        ws = self.treeRoot.find('objective').findall('weight')
        for w in ws:
            weights.append(float(w.get('value')))
        return weights

    def calValueUtils(self):
        """線形効用空間の場合に、各ValueのUtilityを計算する"""
        utils = []
        weights = self.getIssueWeight()

        issues = self.treeRoot.find('objective').findall('issue')
        for i, issue in enumerate(issues):
            tempUtils = [] # issueのValueのUtilities
            maxEval = 1
            values = issue.findall('item')
            for value in values:
                maxEval = max(maxEval, int(value.get('evaluation')))

            for value in values:
                tempUtils.append(weights[i] * float(value.get('evaluation')) / float(maxEval))
            utils.append(tempUtils)
        return utils

    def calBidUtil(self, bid):
        """bidのUtilityを計算"""
        util = 0.0
        issues = self.treeRoot.find('objective').findall('issue')
        for i, issue in enumerate(issues):
            values = issue.findall('item')
            for j, value in enumerate(values):
                if str(bid[i]) == str(value.get('value')):
                    util += self.utils[i][j]
                    break
        return util

    def getValueUtils(self):
        """Value別のUtilityリストを返す"""
        return self.utils

    def getValueUtil(self, issueID, valueID):
        """特定issueの特定valueのUtilityを返す"""
        return self.utils[issueID-1][valueID-1]

    def calAllBidSize(self):
        """全合意案候補数を計算する"""
        size = 1.0
        issues = self.treeRoot.find('objective').findall('issue')
        for issue in issues:
            size *= len(issue.findall('item'))
        return size

    def getAllBidSize(self):
        """全合意案候補数を返す"""
        return self.allBidSize

    def getAllBid(self):
        """全合意案候補を返す"""
        bids = []
        issues = self.treeRoot.find('objective').findall('issue')
        for i, issue in enumerate(issues):
            values = issue.findall('item')
            tempBids = []
            for value in values:
                if len(bids) != 0:
                    for bid in bids:
                        tbid = copy.deepcopy(bid)
                        tbid.append(str(value.get('value')))
                        tempBids.append(tbid)
                else:
                    tempBids.append([str(value.get('value'))])
            bids = tempBids
        return bids

# Preference関連の関数
DEBUG = True

def calDistValueUtil(pref1, pref2, issueID, valueID):
    """特定Valueの効用値の差を計算する"""
    #print issueID, valueID
    return pref1.getValueUtil(issueID, valueID) - pref2.getValueUtil(issueID, valueID)


# 2乗差に基づく効用関数の対立度
def calDistValuesUtil_sq(pref1, pref2, issueID):
    """特定Issueの全Valueの効用値の差を計算する"""
    valueSize1 = pref1.getValueSize(issueID) # 対象となっているIssueのValue数を取得
    if valueSize1 != pref2.getValueSize(issueID):
        sys.stderr.write('2つの効用情報のValueの次元が異なります')
        return None

    valueSize2 = 0 # 対象となっている次のIssueのValue数を取得　(0なら次のIssueは存在しない)
    if issueID != pref1.getIssueSize():
        valueSize2 = pref1.getValueSize(issueID+1)
        if valueSize2 != pref2.getValueSize(issueID+1):
            sys.stderr.write('2つの効用情報のValueの次元が異なります')
            return None

    dist = 0.0
    for valueID1 in range(1, valueSize1+1):
        dist1 = calDistValueUtil(pref1, pref2, issueID, valueID1)
        dist += dist1**2 # 各Valueの差の2乗項
        if valueSize2 != 0:
            for valueID2 in range(1, valueSize2+1):
                dist2 = calDistValueUtil(pref1, pref2, issueID+1, valueID2)
                dist += 2 * dist1 * dist2 / pref1.getValueSize(issueID+1)
    #print dist
    return dist

def calDistPref_sqECO(pref1, pref2):
    """全合意案候補の効用値を合意案候補数で正規化した値を返す"""
    issueSize = pref1.getIssueSize()
    # エラー処理
    if issueSize != pref2.getIssueSize():
        sys.stderr.write('2つの効用情報のIssueの次元が異なります')
        return None

    dist = 0.0
    for issueID in range(1, issueSize+1):
        dist += calDistValuesUtil_sq(pref1, pref2, issueID) / pref1.getValueSize(issueID)
    return dist

def calMOL_sq(pref1, pref2, bids):
    """計算量を無視して、bidsの二乗差の和を返す"""
    if len(bids) == 0: return 0.0 # 対象bidが空集合の場合は0を返す

    dist = 0.0

    # 各bidに関して２乗差を計算し，和をとる
    for bid in bids:
        u  = [pref1.calBidUtil(bid), pref2.calBidUtil(bid)]
        du = u[0] - u[1]
        dist += (pref1.calBidUtil(bid) - pref2.calBidUtil(bid))**2

        if DEBUG:
            print " [DEBUG] ",
            print "(" + str(round(u[0],4)) + " - " + str(round(u[1],4)) + ")^2",
            print " = (" + str(round(du,4)) + ")^2 = " + str(round(du**2,4)),
            print " NOW:" + str(dist)

    ans = dist / len(bids)
    if DEBUG:
        print " [DEBUG] ",
        print "Result: " + str(dist) + " / " + str(len(bids)) + " = " + str(ans)
    return ans

def calDistPref_sq(pref1, pref2):
    """計算量を無視して、単純に全bidの二乗差の和を返す"""
    return calMOL_sq(pref1, pref2, pref1.getAllBid()) # 全bidを取得


def calDistMultiPref_sq(prefs):
    """多者間交渉における全体の効用情報の対立度を返す"""
    dist = 0.0
    bids = prefs[0].getAllBid()
    for bid in bids:
        sumUtil = 0.0
        print "#"
        for pref in prefs:
            print pref.calBidUtil(bid),
            sumUtil += pref.calBidUtil(bid)
            print "(" + str(sumUtil) + ")",
        for pref in prefs:
            dist += (pref.calBidUtil(bid) - sumUtil / len(prefs))**2
    return dist / len(bids)


def getAllPretoBids(pref1, pref2):
    """pref1とpref2の場合のパレート最適Bidのリストを返す"""
    bids1 = pref1.getAllBid()
    bids2 = pref2.getAllBid()

    paretoBids = []
    if bids1 == bids2:
        tempBids = []
        for bid in bids1:
            tempBids.append((pref1.calBidUtil(bid), pref2.calBidUtil(bid), bid))

        # pref1
        tempBids = sorted(tempBids, key=lambda x: float(-x[0])) # 降順ソート
        #print "#1" + str(tempBids)

        maxIdx = -1
        maxTemp = -1.0 # 現時点での相手の効用の最大
        for i, tempBid in enumerate(tempBids):
            # 暫定最大値が更新された場合には、過去に参照した入札、つまり、本エージェントが現時点より高い効用値を望む場合、
            # 全て相手エージェントの効用値が下がることを意味する
            if maxTemp < tempBid[1]:
                # 2つ目のbid以降で，本エージェントが一つ前より効用値が下がる場合
                if i != 0 and tempBids[i-1][0] > tempBid[0]:
                    paretoBids.append(tempBids[maxIdx][2])
                maxTemp = tempBid[1]    # 相手エージェントの暫定最大効用値の更新
                maxIdx  = i             # 相手エージェントの暫定最大効用値をとる添え字の更新

        # pref2
        tempBids = sorted(tempBids, key=lambda x: float(-x[1])) # 降順ソート
        #print "#2" + str(tempBids)

        maxIdx = -1
        maxTemp = -1.0 # 現時点での相手の効用の最大
        for i, tempBid in enumerate(tempBids):
            if DEBUG and i != 0:
                print str(i) + " maxIdx:" + str(maxIdx) + " maxTemp:" + str(maxTemp)
            if maxTemp < tempBid[0]:

                if i != 0 and tempBids[i-1][1] > tempBid[1]:
                    if tempBids[maxIdx][2] not in paretoBids:

                        paretoBids.append(tempBids[maxIdx][2])
                maxTemp = tempBid[0]    # 相手エージェントの暫定最大効用値の更新
                maxIdx  = i             # 相手エージェントの暫定最大効用値をとる添え字の更新

    #print "Pareto: " + str(paretoBids)
    return paretoBids

def calParetoDistPref_sq(pref1, pref2):
    pBids = getAllPretoBids(pref1, pref2)
    dist = 0.0

    for pbid in pBids:
        dist += (pref1.calBidUtil(pbid) -  pref2.calBidUtil(pbid))**2
    return dist / float(len(pBids))


def print3partyDist(prefs):
    """3者の対立度を表示"""

    print "<START> 3-party negotiation: Agent A, Agent B, and Agent C"

    print "Target: All Bids"

    print "[AxB]" + str(calDistPref_sq(prefs[0],prefs[1])),
    print " (ECO Result: " + str(calDistPref_sqECO(prefs[0],prefs[1])) + ")"
    print "[BxC]" + str(calDistPref_sq(prefs[1],prefs[2])),
    print " (ECO Result: " + str(calDistPref_sqECO(prefs[1],prefs[2])) + ")"
    print "[CxA]" + str(calDistPref_sq(prefs[2],prefs[0])),
    print " (ECO Result: " + str(calDistPref_sqECO(prefs[2],prefs[0])) + ")"

    print "[1x2](SQ)" + str(calDistMultiPref_sq( [prefs[0],prefs[1]] ))
    print "[2x3](SQ)" + str(calDistMultiPref_sq( [prefs[1],prefs[2]] ))
    print "[3x1](SQ)" + str(calDistMultiPref_sq( [prefs[2],prefs[0]] ))


    print "Pareto Bids"
    print "[AxB] " + str(calMOL_sq(prefs[0],prefs[1], getAllPretoBids(prefs[0],prefs[1])))
    print "[BxC] " + str(calMOL_sq(prefs[1],prefs[2], getAllPretoBids(prefs[1],prefs[2])))
    print "[CxA] " + str(calMOL_sq(prefs[2],prefs[0], getAllPretoBids(prefs[2],prefs[0])))


# 対象となるxmlのpathをそれぞれ記入 (現時点では3者間交渉にのみ対応)
xmls_path = [
    "./preference/Domain1/Domain1.xml",
    "./preference/Domain2/Domain2_util4.xml",
    "./preference/Domain2/Domain2_util5.xml",
    "./preference/Domain2/Domain2_util6.xml"
]
# それぞれのxmlを読み込み，Preferenceクラスとして格納
prefs = [
    Preference(xmls_path[1]), Preference(xmls_path[2]), Preference(xmls_path[3])
]

DEBUG = False



# ３者間交渉として計算開始
print3partyDist(prefs)

"""
for i in range(3):
    bids = getAllPretoBids(prefs[i], prefs[(i+1)%3])
    print len(bids)
    for j in range(len(bids)):
        print bids[j], prefs[i].calBidUtil(bids[j]), prefs[(i+1)%3].calBidUtil(bids[j])
"""


# TODO: 計算量を考慮した対立度の計算がくそ。デバッグ必要あり
