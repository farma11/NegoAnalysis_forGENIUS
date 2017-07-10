# coding: UTF-8
import re

class Bid(object):
    def __init__(self):
        self.issueSize = 0
        self.valueSize = []




def divideValue(line):
    """文字列Bid[a: a1, b: b2, ...]からValueのListに変換"""
    ansValues = []

    r = re.compile("Bid\[(.*), \]")
    bid = r.search(line)
    if bid != None:
        values = bid.group(1).split(',')
        for value in values:
            r = re.compile(": (.+)")
            v = r.search(str(value))
            if v != None:
                ansValues.append(v.group(1))
    return ansValues
