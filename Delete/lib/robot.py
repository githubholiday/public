# -*- coding:utf-8 -*-
import os
import sys
import requests

__author__ = "Gaoshichen"
__email__ = "shichengao@genome.cn"
__date__ = "2021-12-16"
__version__ = '1.0'
__update__  = '机器人模板'
__doc__ = """
    输入参数：
    title,          str,        标题
    lable_list,     list,       标签，即每列的列名
    content_list,   list,       内容，格式为：[[...],[...],[...]...]
"""


# 义乌小群
#yw_url = "https://qyapi.weixin.qq.com/cgi-bin/webhook/send?key=808541e5-b701-4986-99a5-3439506b480e"
# 分析部大群
#real_url = "https://qyapi.weixin.qq.com/cgi-bin/webhook/send?key=90e20ec2-1f14-434d-83ad-56a0d40e6830"
##纯测序与微生物项目部群
#real_url = "https://qyapi.weixin.qq.com/cgi-bin/webhook/send?key=a5039603-857e-472f-99d0-7022f0a57f10"
real_url = "https://qyapi.weixin.qq.com/cgi-bin/webhook/send?key=a5039603-857e-472f-99d0-7022f0a57f10"
'''
调用方式：
    Robor(title, label_list, content_list, url )
    title:为提醒时的主题
    lable_list：为表头信息，展示时按照|分割展示
    content_list：内容信息，展示时按照|分割展示

'''
class Robots():

    def __init__(self,title,label_list,content_list, url):                                                              #调用该类时，需要传入这3个参数，以下为调用该函数时就会发生的，生成属性和执行方法，而不是用main函数来执行
        self.title = title
        self.label_list = label_list
        self.content_list = content_list
        self.url = url
        self.mes = self.content_list_to_markdown(self.content_list)
        self.label = self.label_list_to_markdown(self.label_list)                                               #to_markdown是后面定义的方法，对输入数据进一步处理
        #self.send_txt(title=self.title,label=self.label,mes=self.mes)
        #self.send_txt(title=self.title,label=self.label,mes=self.mes,url=yw_url)
        self.send_txt(title=self.title,label=self.label,mes=self.mes,url=self.url)
    
    def content_list_to_markdown(self,rets):
        mes=''
        if len(rets) == 0:                                                                                      #若果输入的rets也就时列表是空的，mes不变还是---
            mes = '\t'
        else:
            for ret in rets:
                for i in ret:
                    mes += str(i) + ' | '                                                                     #否则就是把rets列表里的内容用---连接起来
                mes += '\n'                                                           #多的这个'\n|\n'，是因为可能需要Markdown的内容有很多，也就时有多条项目信息要统计，需要换行
        return mes
    
    def label_list_to_markdown(self,rets):
        mes=''
        if len(rets) == 0:
            mes = '\t'
        else:
            for ret in rets:
                mes += str(ret) + ' | '
        return mes

        

    # 机器人提醒,传入参数为str
    def send_txt(self,title,label,mes,url):
        headers = {"Content-Type": "text/plain"}                                            #字典
        send_url = url
        #"""
        if mes == "\t":
            send_data={
                "msgtype": "markdown", #消息类型，固定为markdown
                "markdown":{
                    "content":"# <font color=\"warning\">Sci.db.running Attention!!!!!</font>\n" + 
                    "> <font color=\"info\">**无**</font>"
                }
            }
        else:
            send_data={
                "msgtype": "markdown", # 消息类型，固定为markdown
                "markdown":{
                    "content":"# <font color=\"warning\">Everyone Attention!!!!!</font>\n" +
                    "# <font color=\"warning\">Attention! {0}~ </font>\n".format(title) + "## **{0}**\n{1}".format(label,mes) +
                    "> <font color=\"warning\"> **请注意\n@所有人**</font>"
                }
            }
        res = requests.post(url = send_url, headers = headers, json = send_data)
        
