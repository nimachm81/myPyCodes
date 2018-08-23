## Finance.py


__all__ = ["DatetimeLocalToGMT", "CompanyNews", "Stocks", "StockScreener", "StockReports"]

import numpy as np
import pandas as pd
import datetime
import time
import os
from pandas_datareader import data, wb
import sys
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial
from io import StringIO

import urllib
import re
from bs4 import BeautifulSoup
import requests

import feedparser

from IPython.core.display import display, DisplayObject, HTML

from enum import Enum

from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait # available since 2.4.0
from selenium.webdriver.support import expected_conditions as EC # available since 2.26.0
from selenium.webdriver.common.keys import Keys

import subprocess

from dateutil import tz

from pandas.tseries.offsets import BDay    

import traceback


def DatetimeLocalToGMT(datetime_local):
    zone_gmt = tz.tzutc()
    zone_local = tz.tzlocal()

    time_local = datetime_local.replace(tzinfo=zone_local)

    # Convert time zone
    time_gmt = time_local.astimezone(zone_gmt)
    #time_gmt = time_gmt.replace(tzinfo=zone_local)

    return time_gmt

##----------------------------------------------------------------------------

mp3alarm_path = {'StoreDoor':os.path.join('alarms', 'Store_Door_Chime.mp3'),\
                      'Beep':os.path.join('alarms', 'Beep.mp3'),\
                      'Bleep':os.path.join('alarms', 'Bleep.mp3'),\
                      'Buzz':os.path.join('alarms', 'Buzz.mp3'),\
                      'DoorBuzzer':os.path.join('alarms', 'Door_Buzzer.mp3')}

def SoundAlarm(n=1, p=1, alarm='StoreDoor'):
    for i in range(n):
        os.system('mpg123 {}'.format(mp3alarm_path[alarm]))
        time.sleep(p)
        
##-----------------------------------------------------------------------------

class CompanyNews:
    def __init__(self, update_interval=30):
        self.update_interval = update_interval


    def getReutersCompanyNews(self):
        reuters_rss_url = "http://feeds.reuters.com/reuters/companyNews?format=xml"
            
        feed = feedparser.parse( reuters_rss_url )

        print('feed[ "bozo" ]: ', feed[ "bozo" ])
        print('feed[ "url" ]: ', feed[ "url" ])
        print('feed[ "version" ]: ', feed[ "version" ])
        print('feed[ "channel" ]: ', feed[ "channel" ].keys() )
        print('feed[ "items" ]: ', len(feed[ "items" ]))

        #print(feed[ "items" ][0].keys())
        news_titles = []
        while True:
            feed = feedparser.parse( reuters_rss_url )

            for f in feed[ "items" ]:
                f_title = f['title']
                if f_title not in news_titles:
                    news_titles.append(f_title)
                    print(len(news_titles), " : ", f_title)
                    
                    if f_title.startswith('BRIEF'):
                        SoundAlarm()
                    
            time.sleep(self.update_interval)
        



##-----------------------------------------------------------------------------

class Stocks:
    
    def __init__(self):
        self.stockPrefix = 'stock_'
        self.startDate = datetime.date(2000, 1, 1)

        self.mp3alarm_path = os.path.join('other', 'Store_Door_Chime.mp3')
        return
        
    def GetQuote(self, symbols, site='yahoo'):
        if site=='yahoo':
            N = len(symbols)
            quotes = [None]*N
            for i, name in enumerate(symbols):
                htmlfile = urllib.request.urlopen('http://finance.yahoo.com/quote/{}?p={}'.format(name, name))
                htmltext = htmlfile.read().decode('utf-8')
                regex = ('"ask":{"raw":(.+?),"fmt"')
                pattern = re.compile(regex)
                price = re.findall(pattern, htmltext)
                if isinstance(price, list):
                    quotes[i] = price[0]
                else:
                    quotes[i] = np.nan
            return quotes
            
            
    def GetQuoteBS(self, symbols, site='yahoo'):
        ## TODO: not working.. find the proper scan section in html 
        #BeautifulSoup version 
        if site=='yahoo':
            N = len(symbols)
            quotes = [None]*N
            for i, name in enumerate(symbols):
                url = 'http://finance.yahoo.com/q?s={}'.format(name)
                r = requests.get(url)
                soup = BeautifulSoup(r.text)
                #data = soup.find('span', attrs={'id':'yfs_l84_%s' % name})
                data = soup.find('span', attrs={'id':'yfs_l84_{0}'.format(name)})
                quotes[i] = data.text
            return quotes


    def GetZacksCurrentTimeStamp(self):
        ## This method downloads javascript codes in the html file as well
        from selenium import webdriver

        driver = webdriver.Firefox()
        driver.get('https://www.zacks.com/research/earnings/earnings_display.php')
        htmlpage = driver.page_source


        bs = BeautifulSoup(htmlpage, 'lxml')
        # <input type="hidden" class="hidden" id="current_date" name="current_date" value="1480435725">
        input_timestamp = bs.find(lambda tag: tag.name=='input' \
            and tag.has_attr('id') and tag['id']=='current_date' \
            and tag.has_attr('name') and tag['name']=="current_date"\
            and tag.has_attr('value'))

        if input_timestamp!=None:
            input_timestamp = input_timestamp['value']
        print(input_timestamp)
        driver.close()
        
        self.zackstimestamp = input_timestamp
        return

    def GetZacksLatestSECReports(self, save_file=False, filepath='other/zacksReports.txt'):
        ## https://www.zacks.com/research/earnings/earning_export.php?timestamp=1480441267&tab_id=export_excel   
        user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/53.0.2785.143 Chrome/53.0.2785.143 Safari/537.36'
        headers={'User-Agent':user_agent,} 
        
        html_add = 'https://www.zacks.com/research/earnings/earning_export.php'
        urlrequest = urllib.request.Request(html_add, None, headers)
        htmlfile = urllib.request.urlopen(urlrequest)
        file_txt = htmlfile.read().decode('utf-8')
        
        if save_file:
            with open(filepath, 'w') as file_:
                file_.write(file_txt)
        return file_txt


    def GetZacksLatestSECReports_2(self, save_file=False, filepath='other/zacksReports.txt'):
        ## https://www.zacks.com/research/earnings/earning_export.php?timestamp=1480441267&tab_id=export_excel   
        user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/53.0.2785.143 Chrome/53.0.2785.143 Safari/537.36'
        headers={'User-Agent':user_agent,} 

        html_add = 'https://www.zacks.com/research/earnings/earnings_display.php'
        urlrequest = urllib.request.Request(html_add, None, headers)
        htmlfile = urllib.request.urlopen(urlrequest)
        htmltext = htmlfile.read().decode('windows-1252')
        
        bs = BeautifulSoup(htmltext, 'lxml')
        # <input type="hidden" class="hidden" id="current_date" name="current_date" value="1480435725">
        input_timestamp = bs.find(lambda tag: tag.name=='input' \
            and tag.has_attr('id') and tag['id']=='current_date' \
            and tag.has_attr('name') and tag['name']=="current_date"\
            and tag.has_attr('value'))

        if input_timestamp!=None:
            input_timestamp = input_timestamp['value']
        #print('time stamp :', input_timestamp)        
         
        """    
        regex = ('<input type="hidden" class="hidden" id="current_date" name="current_date" value="(.+?)">')
        pattern = re.compile(regex)
        input_timestamp = re.findall(pattern, htmltext)
        print('time stamp: ', input_timestamp)        
        """
        
        html_add = 'https://www.zacks.com/research/earnings/earning_export.php?timestamp={}&tab_id=export_excel'.format(input_timestamp)
        urlrequest = urllib.request.Request(html_add, None, headers)
        htmlfile = urllib.request.urlopen(urlrequest)
        file_txt = htmlfile.read().decode('utf-8')
        
        if save_file:
            with open(filepath, 'w') as file_:
                file_.write(file_txt)
        return file_txt
        

    def GetZacksLatestSECReportsDataframe(self):
        reports = self.GetZacksLatestSECReports()

        if reports!=None:
            pd_reports = pd.read_csv(StringIO(reports), sep='\t')
            pd_reports = pd_reports.sort_values(by="Report Time", ascending=False)
        
            return pd_reports
        else:
            return None


    def CheckZacksLatestSECReportsRepeatedly(self, intervals_sec=60):
        latest_time = None

        while True:
            reports = self.GetZacksLatestSECReportsDataframe()

            for index, row in reports.iterrows():
                rep_time = datetime.datetime.strptime(row['Report Time'], '%H:%M').time()
                if latest_time==None:
                    latest_time = reports['Report Time'].max()
                    latest_time = datetime.datetime.strptime(latest_time, '%H:%M').time()
                    print('{} :     {},    {} \n{}\n'.format(row['Symbol'], row['Report Time'], row['Company'], row))
                    #self.SoundAlarm()
                else:
                    if rep_time>latest_time:
                        print('{} :     {},    {} \n{}\n'.format(row['Symbol'], row['Report Time'], row['Company'], row))
                        self.SoundAlarm()
                    
            time.sleep(intervals_sec)




    def SoundAlarm(self):
        for i in range(3):
            os.system('mpg123 {}'.format(self.mp3alarm_path))
            time.sleep(1)
                        
            
    def SetCompanyListPaths(self, companyListPath):
        self.companyListPath = companyListPath            

    def SetCompanyDataPath(self, companyDataPath):
        self.companyDataPath = companyDataPath
            
    def GetComapnyList(self, trim=True, sort=True, hasSavedData=False):
        company_list = pd.read_csv(self.companyListPath)
        if trim:
            company_list = company_list.loc[:,['Symbol','Name', 'MarketCap']]
        if sort:
            company_list = company_list.sort_values(by="MarketCap", ascending=False)
        if hasSavedData:
            symbols = company_list['Symbol']
            for i in range(len(symbols)):
                if not self.HasSavedData(symbols[i]):
                    company_list = company_list[company_list.Symbol != symbols[i]]
        return company_list
        
    def HasSavedData(self, symb):
        file_name = self.stockPrefix + symb
        file_path = os.path.join(self.companyDataPath, file_name)
        if os.path.exists(file_path):
            return True
        else:
            return False
            
    def UpdateStockData(self):
        company_list = self.GetComapnyList()
        assert 'Symbol' in company_list.keys()
        
        company_symbols = company_list['Symbol']
        for s in company_symbols:
            symb = s.strip(' ')
            #print(symb, end=' ')
            file_name = self.stockPrefix+ str(symb)

            file_path  = os.path.join(self.companyDataPath, file_name)
            
            if not os.path.exists(file_path):
                endDate = datetime.date.today()
                
                try:
                    df = data.DataReader(symb, 'yahoo', self.startDate, endDate)
                    
                    df.to_hdf(file_path, symb, format='table', mode='w')
                except Exception as e:
                    print(symb, ' error: ', str(e))
                except:
                    print(symb, ' error: ', sys.exc_info()[0])
            else:
                df = pd.read_hdf(file_path, symb)

                date_last = df.index[-1]

                startDate = date_last +  datetime.timedelta(days=1)
                endDate = datetime.date.today()
                diff = 0
                if endDate.weekday()==5:
                    diff = 1
                elif endDate.weekday()==6:
                    diff = 2
                    
                if (startDate +  datetime.timedelta(days=diff)).date() < endDate:
                    try:
                        df = data.DataReader(symb, 'yahoo', startDate, endDate)
                        
                        df.to_hdf(file_path, symb, format='table', mode='a', append=True)
                    except Exception as e:
                        print(symb, ' error: ', str(e))
                    except:
                        print(symb, ' error: ', sys.exc_info()[0])
                        
                
    def plotSymbol(self, symb, duration=[1, 'Y'], curves=['Close', 'Volume'], window=None, figsize=(10, 6),
                    extplt=False, extpltData=[100, 1]):
        """ extpltData = [periode (days), number of extrapolations]
        """
        file_name = self.stockPrefix+ str(symb)
        file_path  = os.path.join(self.companyDataPath, file_name)
        if os.path.exists(file_path):
            df = pd.read_hdf(file_path, symb)
            if window!=None:
                df = df.rolling(window=window).mean()
            
            if duration[1]=='Y':
                len_days = int(duration[0]*365)
                if len_days<len(df):
                    df = df[-len_days:]

            n_p = len(curves)
            f, axarr = plt.subplots(n_p, sharex=True, figsize=figsize)

            p_i = 0
            if 'Close' in curves:
                axarr[p_i].plot(df['Close'])
                if extplt:
                    days_extplt, n_extplt = extpltData
                    for i in range(n_extplt):
                        y = df['Close']
                        y = y[np.logical_not(np.isnan(y))]
                        n_y = (i+1)*days_extplt
                        if n_y<len(y):
                            y_ = y[-n_y:]
                            x_ = np.arange(len(y_))
                            y_L = polynomial.Polynomial(polynomial.polyfit(x_, y_, 1))(x_) 
                            axarr[p_i].plot(df.index[-len(y_):], y_L)
                p_i += 1
            if 'Volume' in curves:
                axarr[p_i].plot(df['Volume']/1.e6)
                p_i += 1
            plt.show()
        else:
            print('Market data for {} not available.'.format(symb))
            


    
    def GetGrowthLinearFit(self, symb, period=30, duration=[1, 'Y']):
        file_name = self.stockPrefix+ str(symb)
        file_path  = os.path.join(self.companyDataPath, file_name)
        if os.path.exists(file_path):
            df = pd.read_hdf(file_path, symb)
            
            slopes = []
            if duration[1]=='Y':
                len_days = int(duration[0]*365)
                if len_days<len(df):
                    df = df[-len_days:]
                    n_pts = int(len_days/period)
                    for i in range(n_pts):
                        y = df['Close']
                        y = y[np.logical_not(np.isnan(y))]
                        n_y = (i+1)*period
                        if n_y<len(y):
                            y_ = y[-n_y:]
                            x_ = np.arange(len(y_))
                            y_L = polynomial.polyfit(x_, y_, 1)
                            assert len(y_L)==2
                            slopes.append(y_L[1])
            return np.array(list(reversed(slopes)))
        else:
            print('Market data for {} not available.'.format(symb))
            return None
            
        
    def GetWeightedAveragedGrowthLinearFit(self, symbs, period=30, duration=[1, 'Y'], weightType='uniform'):
        """ weightType = 
                'uniform' : w = 1
                'linear'  : w = np.arange
        """
        assert weightType in ['uniform', 'linear']
        N = len(symbs)
        res = [None]*N
        for i in range(N):
            symb = symbs[i]
            slopes = self.GetGrowthLinearFit(symb, period, duration)
            assert not np.any(np.isnan(slopes))
            if slopes!=None and len(slopes)>0:
                weight = np.ones(len(slopes))
                if weightType == 'linear':
                    weight = np.arange(len(slopes))+1.0
                weight /= np.sum(weight)/len(slopes)
                avg_weighted = (slopes*weight).mean()
                res[i] = avg_weighted
                assert not np.isnan(res[i])
            else:
                res[i] = np.nan
        return np.array(res)
        
        
        
##--------------------------------------------------------------------------------        
from multiprocessing import Process, Pipe, Queue

class StockScreenerCommands(Enum):
    MonitorYahooGainersLosers = 0


class StockScreener(Process):
    def __init__(self, watchlist=[], commPipe=None, vbose=True):
        
        super(StockScreener, self).__init__()
        
        self.stocksWatchlist = watchlist
        self.watchlistPagePrams = None
        self.YahooGLWatchPrams = None
        self.logFolder = 'Finance'
        
        home = os.path.expanduser('~')
        companylist_path = os.path.join(home, 'Documents', 'MarketData', 'CompanyNames')
        self.NASDAQ_companylist_path = os.path.join(companylist_path, 'NASDAQ_companylist.csv')
        self.NYSE_companylist_path = os.path.join(companylist_path, 'NYSE_companylist.csv')
        self.NYSEMKT_companylist_path = os.path.join(companylist_path, 'AMEX_companylist.csv')
        self.companyDataIntraday_path = os.path.join(home, 'Documents', 'MarketData', 'StocksIntraday') ##where to save
        self.companyDataIntraday30min_path = os.path.join(home, 'Documents', 'MarketData', 'StocksIntraday30min') ##where to save
        
        
        self.df_NASDAQ = pd.read_csv(self.NASDAQ_companylist_path)
        self.df_NYSE = pd.read_csv(self.NYSE_companylist_path)
        self.df_NYSEMKT = pd.read_csv(self.NYSEMKT_companylist_path)

        self.commPipe = commPipe
        
        self.vbose = vbose
        
        self.mp3alarm_path = {'StoreDoor':os.path.join('alarms', 'Store_Door_Chime.mp3'),\
                              'Beep':os.path.join('alarms', 'Beep.mp3'),\
                              'Bleep':os.path.join('alarms', 'Bleep.mp3'),\
                              'Buzz':os.path.join('alarms', 'Buzz.mp3'),\
                              'DoorBuzzer':os.path.join('alarms', 'Door_Buzzer.mp3')}
        
        self.TrackLists = {}

        class paramTypes(Enum):
            intradaylive = 0
            intraday_prev = 1

        self.paramTypes = paramTypes
        
        return
        
    def run(self):
        while True:
            msg = self.commPipe.recv()
            if msg[0] == 'monitorYGL':
                #self.commPipe.send(['resp', 'Monitoring YGL.'])
                gl_symbols = self.GetYahooGainersLosers()
                self.commPipe.send(['data', 'GLSYMBOLS', gl_symbols])
            elif msg[0] == 'quit':
                self.commPipe.send(['resp', 'Quiting process.'])
                return
            else:
                print('Bad command.')
        return
    
    def SoundAlarm(self, n=1, p=1, alarm='StoreDoor'):
        for i in range(n):
            os.system('mpg123 {}'.format(self.mp3alarm_path[alarm]))
            time.sleep(p)
    
        
    def CloseDriver(self):
        if self.watchlistPagePrams!=None:
            driver = self.watchlistPagePrams['driver']
            driver.quit()
        
        
    def AddToWatchlist(self, watchlist):
        for w in watchlist:
            if w not in self.stocksWatchlist:
                self.stocksWatchlist.append(w)

    def StartGoogleFinanceWatchlist(self):
        # Create a new instance of the Chrome/Firefox driver
        driver = webdriver.Chrome()
        #driver = webdriver.Firefox()
        #driver = webdriver.Remote("http://localhost:4444/wd/hub", webdriver.DesiredCapabilities.HTMLUNIT.copy())
        #driver = webdriver.Remote("http://localhost:4444/wd/hub", webdriver.DesiredCapabilities.HTMLUNITWITHJS)
        #driver = webdriver.Remote('http://127.0.0.1:4444/wd/hub', webdriver.DesiredCapabilities.CHROME)

        stocksWatchlist = self.stocksWatchlist

        n_stocks  = len(stocksWatchlist)
        stocksPrices = [None]*n_stocks
        windowNames = [None]*n_stocks
        windowTitles = [None]*n_stocks
        pageLoaded = [None]*n_stocks
        win_id = None

        # go to the google home page
        #driver.set_window_size(1920,1080)

        n_load_try = 3
        for i in range(n_stocks):
            try:
                if i>0:
                    last_handle = driver.window_handles[-1]
                    while driver.window_handles[-1]==last_handle:
                        p = subprocess.Popen(['xdotool', 'windowfocus', win_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        out, err = p.communicate()
                        #print(out.decode('utf-8'))
                    
                        driver.find_element_by_tag_name('body').send_keys(Keys.CONTROL + 't') 
                        #driver.find_element_by_tag_name('body').send_keys(Keys.chord(Keys.CONTROL, 't'))
                        
                        #time.sleep(0.5)
                    driver.switch_to.window(driver.window_handles[-1])
                    
                
                if i==0:
                    driver.get(r"about:blank")
                    p = subprocess.Popen(['xdotool', 'getactivewindow'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    win_id = out.decode('utf-8')
                
                driver.get(r"http://www.google.com/finance")

                windowNames[i] = driver.window_handles[-1]
                time.sleep(0.2)
                well_loaded = False
                for j in range(n_load_try):
                    enteredText = None
                    inputElement = driver.find_element_by_name("q")
                    while enteredText!=stocksWatchlist[i]:
                        inputElement.clear()
                        inputElement.send_keys(stocksWatchlist[i])
                        time.sleep(0.2)
                        enteredText = inputElement.get_attribute('value')
                    inputElement.submit()
                    #print('entered field: ', enteredText)

                    #WebDriverWait(driver, 10).until(EC.title_contains(stocksWatchlist[i]))
                    windowTitles[i] = driver.title

                    ##TODO: double check symbol, reload if necessary        
                    company_info = driver.find_element_by_css_selector('div.appbar-center')
                    #print(company_info)
                    company_name = company_info.find_element_by_css_selector('div.appbar-snippet-primary').text.strip()
                    company_symb = company_info.find_element_by_css_selector('div.appbar-snippet-secondary').text.strip()
                    #print('Name: {}, Symbol: {} '.format(company_name, company_symb) )
                    ## TODO: if name or symbol is not correct repeat..

                    if company_symb=='({})'.format(stocksWatchlist[i]):
                        if self.vbose:
                            print('{}: {} loaded successfully!'.format(stocksWatchlist[i], company_name))
                        well_loaded = True
                        break
                    else:
                        if self.vbose:
                            print('{}: {} not loaded!'.format(j, stocksWatchlist[i]))
                pageLoaded[i] = well_loaded
                time.sleep(0.3)
                        
                
            except:
                print('Got an exception: ', sys.exc_info()[0])
                driver.save_screenshot(os.path.join(self.logFolder, 'screenshot.png'))

        #print(windowTitles)
        self.watchlistPagePrams = {'driver':driver, 'win_id':win_id, 'stocksPrices':stocksPrices, \
            'windowNames':windowNames, 'windowTitles':windowTitles}
        
        return
        
                  
    def getFlashChartData(self, driver):
        #chartElement = driver.find_element_by_id("chart_anchor")
        chartElement = driver.find_element_by_id("chart-section")
        chartHtml = chartElement.get_attribute('outerHTML')
        if self.vbose:
            print(chartHtml)
        chart_str = chartHtml #BeautifulSoup(chartHtml, 'lxml').text.strip()
        chart_str = driver.page_source
        chart_st = chart_str.find(r'chart:')
        chart_end = chart_str.find(r'}', chart_st)
        chart = chart_str[chart_st: chart_end+1]
        ##find _5d:" ... "
        data_5d_st = chart.find(r'_5d:"')+5
        data_5d_end = chart.find(r'"', data_5d_st)
        data_5d = chart[data_5d_st: data_5d_end]
        ##find the last \na 
        inds_na = []
        while True:
            ind_na = None
            if len(inds_na)>0:
                ind_na = data_5d.find(r'\na', inds_na[-1]+1)
            else:
                ind_na = data_5d.find(r'\na')
            if ind_na>=0:
                inds_na.append(ind_na)
            else:
                break
                
        data = data_5d[inds_na[-1]:-1]
        return data
              
 
 
    ##symb=NASDAQ:AMD, intval=60, period=25
    def getIntradayOrHistoricPrices(self, symb, intval='60', period='25m', extendex_hours=False):
        ## interval:seconds, period:m: minutes  d: days   Y: years
        assert int(intval)%30==0
        assert period[-1] in ['m', 'd', 'Y']
        
        mkt_ind = symb.find(r':')
        assert mkt_ind >=0
        mkt = symb[:mkt_ind]
        symb = symb[mkt_ind+1:]
        
        if mkt=='NASDAQ':
            mkt = 'NASD'
        url = 'https://finance.google.ca/finance/getprices?q={}&x={}&i={}&p={}'.format(symb, mkt, intval, period)
        
        if extendex_hours:
            url += '&sessions=ext_hours'
        
        html_req = urllib.request.urlopen(url)
        datatext = html_req.read().decode('utf-8')
        #print(url)
        #print(datatext)
        
        assert datatext.find('COLUMNS=DATE,CLOSE,HIGH,LOW,OPEN,VOLUME')>=0
        assert datatext.find('INTERVAL={}'.format(intval))>0
        
        data_st = datatext.find(r'DATA=')
        data_st = datatext.find('\n', data_st)+1
        
        datatext = datatext[data_st:]
        
        ##remove TIME_OFFSET=xx\n
        while True:
            to_0 = datatext.find('TIMEZONE_OFFSET')
            if to_0<0:
                break
            to_1 = datatext.find('\n', to_0)
            assert to_1>=0
            datatext = datatext[:to_0] + datatext[to_1+1:]
            
        ##
        data = 'DATE,CLOSE,HIGH,LOW,OPEN,VOLUME\n'+datatext
        #print(data)
        
        df = pd.read_csv(StringIO(data), sep=',')
        if len(df)==0:
            return df
            
        assert df.loc[0, 'DATE'].startswith('a')

        date_0 = None
        for i in range(len(df)):
            if df.loc[i, 'DATE'].startswith('a'):
                date_0 = int(df.loc[i, 'DATE'].replace('a', ''))
                df.loc[i, 'DATE'] = i
                
            df.loc[i, 'DATE'] = datetime.datetime.fromtimestamp(date_0 + int(df.loc[i, 'DATE'])*int(intval))
        #print(df)
        
        return df

    def getLastClose(self, symb):
        df = self.getIntradayOrHistoricPrices(symb, intval='86400', period='10d')
        if len(df)==0:
            return None
        ind_last = len(df)-1
        return df.loc[ind_last, 'CLOSE']
        
    def getLastPrice(self, symb):
        df = self.getIntradayOrHistoricPrices(symb, intval='30', period='25m')
        if len(df)==0:
            return None
        ind_last = len(df)-1
        return df.loc[ind_last, 'CLOSE']
        
    def getOpenPriceFast(self, symb):
        ##not valid after more than 50mins of market open
        df = self.getIntradayOrHistoricPrices(symb, intval='30', period='50m')
        if len(df)==0:
            return None
        return df.loc[0, 'OPEN']
    
    def getAverageHistoricPriceAndVolume(self, symb, n_day=30):        
        df = self.getIntradayOrHistoricPrices(symb, intval='1800', period='{}d'.format(n_day))
        if len(df)==0:
            return None
        price_avg = df['CLOSE'].mean()
        vol_avg = df['VOLUME'].mean()
        return price_avg, vol_avg

    def updatePrices(self, every_n_sec=5, refresh_page_every_n_min=5):
        time_last = time.perf_counter()
        page_refreshed = []
        
        refresh_page_every_n_sec = refresh_page_every_n_min*60
        
        n_stocks = len(self.stocksWatchlist)
        
        driver = self.watchlistPagePrams['driver']
        stocksPrices = self.watchlistPagePrams['stocksPrices']
        windowNames = self.watchlistPagePrams['windowNames']

        try:
            while True:
                just_refreshed = False
                for i in range(n_stocks):
                    driver.switch_to.window(windowNames[i])
                    priceElement = driver.find_element_by_id("price-panel")
                    price_text = priceElement.text
                    #print(price_text)
                    priceHtml = priceElement.get_attribute('outerHTML')
                    #print(priceHtml)

                    bs = BeautifulSoup(priceHtml, 'lxml')
                    data = bs.find('span', attrs={'class':'pr'}).text.strip()
                    stocksPrices[i] = data
                    
                    ##TODO: update based on the realtime clock on the page
                    tima_elapsed = time.perf_counter() - time_last
                    if tima_elapsed>refresh_page_every_n_sec and not just_refreshed:
                        if i not in page_refreshed:
                            driver.refresh()
                            page_refreshed.append(i)
                            just_refreshed = True
                        if len(page_refreshed)==n_stocks:
                            time_last = time.perf_counter()
                            page_refreshed = []
                            
                if self.vbose:
                    print(stocksPrices)
                time.sleep(every_n_sec)
                just_refreshed = False
                #break
        except Exception as err:
            print('Got an exception: ', str(err), err.args)
            driver.save_screenshot(os.path.join(self.logFolder, 'screenshot.png'))
        except:
            print('Got an exception: ', sys.exc_info()[0])
            driver.save_screenshot(os.path.join(self.logFolder, 'screenshot.png'))
        finally:
            #driver.quit()
            pass


    def FindMarketAndSymbol(self, symb, name):
        df = self.df_NASDAQ
        res = df.loc[df['Symbol'] == symb]
        if len(res)>0 and self.NormalizeCompanyName(res['Name'].iloc[0])==self.NormalizeCompanyName(name):
            return r'NASDAQ:{}'.format(symb)
        df = self.df_NYSE
        res = df.loc[df['Symbol'] == symb]
        if len(res)>0 and self.NormalizeCompanyName(res['Name'].iloc[0])==self.NormalizeCompanyName(name):
            return r'NYSE:{}'.format(symb)
        df = self.df_NYSEMKT
        res = df.loc[df['Symbol'] == symb]
        if len(res)>0 and self.NormalizeCompanyName(res['Name'].iloc[0])==self.NormalizeCompanyName(name):
            return r'NYSEMKT:{}'.format(symb)
        if self.vbose:
            print('{}:{} not found.'.format(symb, name))
        return None

    def NormalizeCompanyName(self, name):
        name_norm = name.lower().replace(',', '').replace('.', '').replace('inc', '')\
            .replace('limited', '').replace('ltd', '').replace('llc', '')\
            .replace('plc', '').replace('corp', '').strip()
        return name_norm


    def ExtractYahooGainerSymbols(self, table):
        #td class="data-col0 Pstart(10px)"
        bs = BeautifulSoup(table, 'lxml')
        col_0 = bs.findAll("td", class_="data-col0" )
        symbols_dic = {}
        #print(col_0)
        
        
        for i in range(len(col_0)):
            symb_i = col_0[i].find('div', class_=r"Ell").text
            name_i = col_0[i].find('div', class_=r"Fz(xs)").text
            symbols_dic[symb_i] = name_i
        #print(symbols)
        #print(names)
        return symbols_dic

    def ExtractYahooGainerDataframe(self, table):
        dfs = pd.read_html(table)
        assert len(dfs)==1
        df = dfs[0]
        
        """
        symbols = [df.loc[i, 'Symbol'] for i in range(len(df))]
        symbols_init = [symbols[i].split()[0] for i in range(len(symbols))]
        symbols_init_up = ["".join([l for l in symbols_init[i] if l.isupper()]) for i in range(len(symbols_init))]
        symbols_init_up = [symbols_init_up[i][0:len(symbols_init_up[i])-1] for i in range(len(symbols_init_up))]
        """
        return df


    def GetYahooGainersOrLosers(self, gainersOrLosers='gainers'):
        assert gainersOrLosers in ['gainers', 'losers']
        driver = webdriver.Chrome()
        
        driver.get(r'about:blank'.format(gainersOrLosers))
        ##----
        p = subprocess.Popen(['xdotool', 'getactivewindow'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        win_active = out.decode('utf-8')
        
        p = subprocess.Popen(['xdotool', 'windowminimize', win_active], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        ##---

        driver.get(r'http://finance.yahoo.com/{}'.format(gainersOrLosers))
        
        if gainersOrLosers=='gainers':
            WebDriverWait(driver, 10).until(EC.title_contains('Gainers'))
        else:
            WebDriverWait(driver, 10).until(EC.title_contains('Losers'))
        tables = driver.find_elements_by_tag_name('table')

        #if 'sorry, but we were unable to retrieve your list' in driver.page_source:
        #    print('Table not updated yet.')

        symbs_list = []
        for table in tables:
            header = ''
            try:
                headers = table.find_elements_by_tag_name('th')
                for h in headers:
                    header += h.text+' '
                header = header.strip()
                #print(header)
            except:
                pass
            if 'Symbol' in header and 'Last Price' in header and 'Market Time' in header:
                symbols = self.ExtractYahooGainerSymbols(table.get_attribute('outerHTML'))

                for symb, name in symbols.items():
                    mkt_symb = self.FindMarketAndSymbol(symb, name)
                    if mkt_symb!=None:
                        symbs_list.append(mkt_symb)

        driver.quit()
        return symbs_list
       

    def SetupYahooGainerLoserWatch(self, onlyGainers=False):
        driver = webdriver.Chrome()

        n_gl = 2
        if onlyGainers:
            n_gl = 1
        windowNames = [None]*n_gl
        windowTitles = [None]*n_gl
        pageLoaded = [None]*n_gl
        win_id = None
        
        gainersOrLosers = ['gainers', 'losers']
        if onlyGainers:
            gainersOrLosers = ['gainers']

        n_load_try = 3
        for i in range(n_gl):
            try:
                if i>0:
                    last_handle = driver.window_handles[-1]
                    while driver.window_handles[-1]==last_handle:
                        p = subprocess.Popen(['xdotool', 'windowfocus', win_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        out, err = p.communicate()
                    
                        driver.find_element_by_tag_name('body').send_keys(Keys.CONTROL + 't') 
                        
                    driver.switch_to.window(driver.window_handles[-1])
                
                if i==0:
                    driver.get(r"about:blank")
                    p = subprocess.Popen(['xdotool', 'getactivewindow'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    win_id = out.decode('utf-8')
                
                driver.get(r'http://finance.yahoo.com/{}'.format(gainersOrLosers[i]))

                windowNames[i] = driver.window_handles[-1]
                windowTitles[i] = driver.title
            except:
                print('Got an exception: ', sys.exc_info()[0])
                driver.save_screenshot(os.path.join(self.logFolder, 'screenshot.png'))

        #print(windowTitles)
        self.YahooGLWatchPrams = {'driver':driver, 'win_id':win_id, \
            'windowNames':windowNames, 'windowTitles':windowTitles}
        
        return
            

    def GetYahooGainersLosers(self, refresh=True, onlyGainers=True):
        if self.YahooGLWatchPrams==None:
            self.SetupYahooGainerLoserWatch(onlyGainers=onlyGainers)
            
        n_gl = 2
        if onlyGainers:
            n_gl = 1
        driver = self.YahooGLWatchPrams['driver']
        windowNames = self.YahooGLWatchPrams['windowNames']
        gl_symbols = [None]*n_gl
        
        gainersOrLosers = ['gainers', 'losers']
        if onlyGainers:
            gainersOrLosers = ['gainers']
        
        in_titles = ['Gainers', 'Losers']
        if onlyGainers:
            in_titles = ['Gainers']

        try:
            for i in range(n_gl):
                driver.switch_to.window(windowNames[i])
                if refresh:
                    driver.refresh()
                
                title = driver.title
                if in_titles[i] not in title:
                    driver.get(r'http://finance.yahoo.com/{}'.format(gainersOrLosers[i]))
                
                WebDriverWait(driver, 5).until(EC.title_contains(in_titles[i]))        
                tables = driver.find_elements_by_tag_name('table')
                
                symbs_list = []
                for table in tables:
                    header = ''
                    try:
                        headers = table.find_elements_by_tag_name('th')
                        for h in headers:
                            header += h.text+' '
                        header = header.strip()
                        #print(header)
                    except:
                        pass
                    if 'Symbol' in header and 'Last Price' in header and 'Market Time' in header:
                        symbols = self.ExtractYahooGainerSymbols(table.get_attribute('outerHTML'))
                        
                        for symb, name in symbols.items():
                            mkt_symb = self.FindMarketAndSymbol(symb, name)
                            if mkt_symb!=None:
                                symbs_list.append(mkt_symb)
                
                gl_symbols[i] = symbs_list
                    
        except Exception as err:
            print('GetYahooGainersLosers: ', str(err), err.args)
            driver.save_screenshot(os.path.join(self.logFolder, 'screenshot.png'))
        except:
            print('GetYahooGainersLosers: ', sys.exc_info()[0])
            driver.save_screenshot(os.path.join(self.logFolder, 'screenshot.png'))
        
        return gl_symbols
            


    def GetNotableGainerLoserStocks(self, symbols, price_range=[1.0, 5.0], volPrice_avg=10000, close_open_jump_perc_max=10):
        ##less than 10$
        ##small premarket change
        ##breaking 20-week average sharply (by 10%) in the last few days
        ##to report soon
        ##real-time data available
        ##decent liquidity
        
        ##TODO: keep if traded volume much higher than average
        
        symbs_note = []
        for symb in symbols:
            last_close =  self.getLastClose(symb) 
            last_price = self.getLastPrice(symb)
            ##average volume/30mins
            price_vol_avg_30d = self.getAverageHistoricPriceAndVolume(symb, n_day=30)
            open_price = self.getOpenPriceFast(symb)
            
            if self.vbose:
                print('{}: Close:{} Open:{} '.format(symb, last_close, open_price), end=' ')

            if last_close==None or last_price==None or price_vol_avg_30d==None or open_price==None:
                continue
            
            price_avg_30d, vol_avg_30d = price_vol_avg_30d
            
            if self.vbose:
                print('PV: ',  vol_avg_30d*price_avg_30d)
            
            premarket_change = (open_price-last_close)/last_close*100
            
            notable = True
            
            if not (last_price>price_range[0] and last_price<price_range[1]):
                notable = False
            if not ( abs(premarket_change)<close_open_jump_perc_max ):
                notable = False
            if not (vol_avg_30d*price_avg_30d)>volPrice_avg:
                notable = False
                
            if notable:
                symbs_note.append(symb)
                if self.vbose:
                    print(' !!!!! ')
                
        return symbs_note

    
    def SetStocksToTrack(self, list_name, stock_list):
        stock_list = list(set(stock_list))
        self.TrackLists[list_name] = {}
        stocksToTrack = stock_list
        stocksToTrackData = [None]*len(stock_list)
        stocksToTrackNextUpdateTime = [None]*len(stock_list)
        stocksToTrackIndicators = [None]*len(stock_list)
        for i in range(len(stock_list)):
            self.stocksToTrackIndicators[i] = {}
        stocksRemovedFromTrackList = []
        
        self.TrackLists[list_name]['stocksToTrack'] = stocksToTrack
        self.TrackLists[list_name]['stocksToTrackData'] = stocksToTrackData
        self.TrackLists[list_name]['stocksToTrackNextUpdateTime'] = stocksToTrackNextUpdateTime
        self.TrackLists[list_name]['stocksToTrackIndicators'] = stocksToTrackIndicators
        self.TrackLists[list_name]['stocksRemovedFromTrackList'] = stocksRemovedFromTrackList
        return
        
    def AddStockToTrackList(self, list_name, stock):
        if isinstance(stock, list):
            for s in stock:
                self.AddStockToTrackList(list_name, s)
            return
        
        assert list_name in self.TrackLists
        trackList = self.TrackLists[list_name]
        stocksToTrack = trackList['stocksToTrack']
        stocksToTrackData = trackList['stocksToTrackData']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        stocksToTrackNextUpdateTime = trackList['stocksToTrackNextUpdateTime']
        stocksToTrackIndicators = trackList['stocksToTrackIndicators']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        
        if stock not in stocksToTrack  \
                and stock not in stocksRemovedFromTrackList:
            stocksToTrack.append(stock)
            stocksToTrackData.append(None)
            stocksToTrackNextUpdateTime.append(None)
            stocksToTrackIndicators.append({})

    def RemoveStockFromTrackList(self, list_name, stock):
        if isinstance(stock, list):
            for s in stock:
                self.RemoveStockFromTrackList(list_name, s)
            return

        assert list_name in self.TrackLists
        trackList = self.TrackLists[list_name]
        stocksToTrack = trackList['stocksToTrack']
        stocksToTrackData = trackList['stocksToTrackData']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        stocksToTrackNextUpdateTime = trackList['stocksToTrackNextUpdateTime']
        stocksToTrackIndicators = trackList['stocksToTrackIndicators']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']

        if stock in stocksToTrack:
            ind = stocksToTrack.index(stock)
            del stocksToTrack[ind]
            del stocksToTrackData[ind]
            del stocksToTrackNextUpdateTime[ind]
            del stocksToTrackIndicators[ind]
            stocksRemovedFromTrackList.append(stock)
            print('{} removed from {}!'.format(stock, list_name))
            

    def UpdateTrackedStockData(self, list_name):
        assert list_name in self.TrackLists
        trackList = self.TrackLists[list_name]
        stocksToTrack = trackList['stocksToTrack']
        stocksToTrackData = trackList['stocksToTrackData']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        stocksToTrackNextUpdateTime = trackList['stocksToTrackNextUpdateTime']
        stocksToTrackIndicators = trackList['stocksToTrackIndicators']

        t_now = datetime.datetime.now()
        for i in range(len(stocksToTrack)):
            t_update = stocksToTrackNextUpdateTime[i]
            if t_update==None or t_update<t_now:
                ##update
                symb =  stocksToTrack[i]
                df =  stocksToTrackData[i]
                period = '25m'
                intval = '60'
                merge_df = True
                if not isinstance(df, pd.DataFrame):
                    merge_df = False
                    period = '1d'
                df_new = None
                new_data_available = False
                try:
                    df_new = self.getIntradayOrHistoricPrices(symb, intval, period)
                    new_data_available = True
                except:
                    pass
                if new_data_available:
                    assert isinstance(df_new, pd.DataFrame)
                    assert 'DATE' in df_new
                    df_new = df_new.set_index('DATE')
                    if merge_df:
                        ##merge
                        index = df_new.index
                        for j in range(len(index)):
                            df.loc[index[j]] = df_new.loc[index[j]]
                        
                    else:
                        t_open = t_now.replace(hour=9, minute=30, second=0, microsecond=0)
                        df = df_new[df_new.index>=t_open]
                
                if isinstance(df, pd.DataFrame):
                    df.sort_index(ascending=True,  inplace=True)
                    stocksToTrackData[i] = df
                    ##schedule next update
                    index = df.index
                    if len(index)>0:
                        t_last = index[-1]
                        dt = datetime.timedelta(minutes=1, seconds=5)
                        if len(index)>1:
                            dt = t_last - index[-2] +  datetime.timedelta(seconds=5)
                        t_next = t_last + dt
                        stocksToTrackNextUpdateTime[i] = t_next
                if self.vbose:
                    print('{} - {}'.format(symb, stocksToTrackNextUpdateTime[i]))
                    print(stocksToTrackData[i])
                    print('-'*60)
                if isinstance(stocksToTrackData[i], pd.DataFrame):
                    print('{} ({}-{})'.format(i, symb, len(stocksToTrackData[i])), end=' ')
                else:
                    print('{} ({}-None)'.format(i, symb), end=' ')
                
        print('')
        if self.vbose:
            print('='*60)
        return
 
 
    ##TODO: correct aftehour mistakes
    def GetLastDayTradingInfo(self, symb):
        dir_path = self.companyDataIntraday30min_path
        if self.HasSavedData(symb, dir_path):
            assert symb.find(':')>=0
            symb_bare = symb[symb.find(':')+1:]

            file_name = symb.replace(':', '_')
            file_path = os.path.join(dir_path, file_name)
          
            store = pd.HDFStore(file_path)
            keys = store.keys()
            assert r'/'+symb_bare in keys

            nrows = store.get_storer(symb_bare).nrows
            
            start = None
            if nrows>1:
                start = nrows-1
            #df = pd.read_hdf(file_path, symb, start=start)
            df = store.select(symb_bare, start=start)
 
            date_last = self.NpPdDatetimeToPyDatetime(df.index[-1]).date()
            df = store.select(symb_bare, where='index>=date_last')
            
            store.close()
            df = df[np.logical_and(df.index.time>=datetime.time(9, 30), df.index.time<=datetime.time(16, 0))]
            
            df = df.sort_index(ascending=True)
            
            trading_period_mins = (self.NpPdDatetimeToPyDatetime(df.index[-1]) - \
                        self.NpPdDatetimeToPyDatetime(df.index[0])).seconds//60
            trading_volume = np.sum(df['VOLUME'])
            trading_price_avg = np.mean(df['CLOSE'])
            trading_price_close = df['CLOSE'][-1]
            
            return {'dataframe':df, 'period':trading_period_mins, \
                'volume':trading_volume, 'price-average':trading_price_avg, \
                'price-close':trading_price_close}
        return None
        
        

    def GetIntraDayMinuteTradingInfo(self, symb, date, mins_max=60):
        """ mins_max: provide 1-mins_max minute average and maximum data
        """
        dir_path = self.companyDataIntraday_path
        if self.HasSavedData(symb, dir_path):
            assert symb.find(':')>=0
            symb_bare = symb[symb.find(':')+1:]

            file_name = symb.replace(':', '_')
            file_path = os.path.join(dir_path, file_name)
          
            store = pd.HDFStore(file_path)
            keys = store.keys()
            assert r'/'+symb_bare in keys

            nrows = store.get_storer(symb_bare).nrows
            
            start = None
            if nrows>1:
                start = nrows-1
            #df = pd.read_hdf(file_path, symb, start=start)
            df = store.select(symb_bare, start=start)
            
            df = self.CorrectMarketDataErrors(df)
 
            day_open = datetime.datetime.combine(date, datetime.time(hour=9, minute=30))
            day_close = day_open.replace(hour=16, minute=0, second=0, microsecond=0)

            df = store.select(symb_bare, where='index>=day_open')
            
            store.close()
            df = df[np.logical_and(df.index>=day_open, df.index<=day_close)]
            
            df = df.sort_index(ascending=True)
            
            P = df['CLOSE'].values
            V = df['VOLUME'].values
            T = df.index.values
            N = len(T)
            if N<10:
                return None
            
            dT = T[1:] - T[:-1]
            dT_mins = self.NpPdtimedeltaArrayToMins(dT)
            
            V_ = np.concatenate((np.zeros(1), V[1:]))
            dT_ = np.concatenate((np.zeros(1), dT_mins))
            
            m_max = 60
            V_max = np.zeros(m_max)
            V_avg = np.zeros(m_max)
            ##TODO: take into account 2min series
            for i in range(1, m_max):
                V_i = np.cumsum(V_)[i:] - np.cumsum(V_)[:-i] 
                dT_i = np.cumsum(dT_)[i:] - np.cumsum(dT_)[:-i]
                #print(dT_i)
                dT_eq_i = (dT_i==i)
                if np.sum(dT_eq_i)>50:
                    V_i_max = (V_i*dT_eq_i).max()
                    V_i_avg = (V_i*dT_eq_i).mean()
                    V_max[i] = V_i_max
                    V_avg[i] = V_i_avg
                
                
            return {'dataframe':df, \
                'volume-avg':V_avg, 'volume-max':V_max, 'price-average':P.mean(), \
                'price-max':P.max(), 'price-min':P.min(), 'price-open':P[0], 'price-close':P[-1]}
        return None




    """ dv>0 dp>0  ---> rush to buy
        dv>0 dp<0  ---> rush to sell in a liquid market
        dv<0 dp<0  ---> rush to sell, market illiquid
        dv<0 dp>0  ---> confusion?
        
        dp>2%, dv>0 in the last n minutes   
    """
    ##TODO: correct afterhour mistakes
    ##TODO: 1-2 min opening volume much higher than last-day average volume and price increasing
    ##TODO: 1-2 min volume should be twice the average 10 minute volume
    ##TODO: >10 min trend should have continued with the same rate during the last hour
    ##TODO: find piecewise linear fit or second order fit and detect the minimums 
    ##TODO: if price increases in moderate volume and then stays flat on low volume, the price normally increases

    ##TODO: when the price is decreasing in increasing volume it will jump back up
    ##TODO: when the price is increasing in increasing volume it will jump back down
    ##TODO: rapid movements after open and befor close each day
    def GetIndicatorsForTrackedStocks(self, list_name, pv_th_1=30.0e3, pv_th_60=500.0e5, dp_th=0.01, d_pv_th=1.0):
        ## dp/dt > 2%
        ## dv/dt > 0
        ## v*t>100*p_budget
        ## dt: 1, 2, 10, 20, 30, 60
        assert list_name in self.TrackLists
        trackList = self.TrackLists[list_name]
        stocksToTrack = trackList['stocksToTrack']
        stocksToTrackData = trackList['stocksToTrackData']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        stocksToTrackNextUpdateTime = trackList['stocksToTrackNextUpdateTime']
        stocksToTrackIndicators = trackList['stocksToTrackIndicators']

        for i in range(len(stocksToTrack)):
            symb =  stocksToTrack[i]
            df =  stocksToTrackData[i]
            if not isinstance(df, pd.DataFrame):
                continue
            if len(df)<=2:
                continue
            N = len(df)
            P = df['CLOSE'].values
            V = df['VOLUME'].values
            T = df.index.values
            
            d_pv = {}
            d_pv_F = {}
            
            t_n = T[-1]
            #print((T[-1]-T[-2]).astype('timedelta64[m]').item())
            for dt_i in range(1, 61):
                #dt = datetime.timedelta(minutes=dt_i)
                dt = np.timedelta64(dt_i, 'm')
                ind_t0, ind_t1 = None, None
                for j in range(N-2, -1, -1):
                    if ind_t1==None and (t_n-T[j])>=dt:
                        ind_t1 = j
                    elif ind_t0==None and ind_t1!=None and (T[ind_t1]-T[j])>=dt:
                        ind_t0 = j
                        break
                if ind_t0==None or ind_t1==None:
                    continue
                
                d_t0t1 = (T[ind_t1]-T[ind_t0]).astype('timedelta64[m]').item().total_seconds()/60
                d_t1tn = (t_n-T[ind_t1]).astype('timedelta64[m]').item().total_seconds()/60
                v_t0t1 = 0.0
                p_t0t1 = 0.0
                for j in range(ind_t0+1, ind_t1+1):
                    v_t0t1 += V[j]*(T[j]-T[j-1]).astype('timedelta64[m]').item().total_seconds()/60
                    p_t0t1 += P[j]*(T[j]-T[j-1]).astype('timedelta64[m]').item().total_seconds()/60
                p_t0t1 /= d_t0t1
                v_t0t1 /= d_t0t1
                
                v_t1tn = 0.0
                p_t1tn = 0.0
                for j in range(ind_t1+1, N):
                    v_t1tn += V[j]*(T[j]-T[j-1]).astype('timedelta64[m]').item().total_seconds()/60
                    p_t1tn += P[j]*(T[j]-T[j-1]).astype('timedelta64[m]').item().total_seconds()/60
                p_t1tn /= d_t1tn
                v_t1tn /= d_t1tn
                
                d_pv[dt_i] = [t_n, (p_t1tn*v_t1tn - p_t0t1*v_t0t1)/(p_t0t1*v_t0t1)]
                
                pv_th = (pv_th_60 - pv_th_1)*(dt_i-1)/59 + pv_th_1
                d_pv_F[dt_i] = [t_n, (p_t1tn*v_t1tn - p_t0t1*v_t0t1)/(p_t0t1*v_t0t1)\
                    *((p_t1tn*v_t1tn)>pv_th)*((p_t1tn - p_t0t1)/p_t0t1>dp_th)]
            
            d_pv_max = 0.0
            dt_max = None
            date_max = None
            for dt_i in d_pv_F:
                if d_pv_F[dt_i]!=None and d_pv_F[dt_i][1]>d_pv_max:
                    d_pv_max = d_pv_F[dt_i][1]
                    date_max = d_pv_F[dt_i][0]
                    dt_max = dt_i
            if d_pv_max>d_pv_th:
                print('-'*60+ '\n{} : time: {} DT: {} indicator: {}\n'.format(symb, date_max, dt_max, d_pv_max) + '-'*60)
            
            if dt_max!=None and dt_max<=2:
                lastDayInfo = self.GetLastDayTradingInfo(symb)
                if lastDayInfo!=None:
                    vol_avg = lastDayInfo['volume']
                    period_mins = lastDayInfo['period']
                    vol_avg_dt = vol_avg/period_mins*dt_max
                    
                    if V[-1]>50.0*vol_avg_dt:
                        print('!'*50)
                        print('-'*60+ '\n{} : time: {} DT: {} indicator: {}\n'.format(symb, date_max, dt_max, d_pv_max) + '-'*60)
                        
            
            stocksToTrackIndicators[i]['d_pv'] = d_pv
        
        return


    
    def ApplyRealtimeModelToTrackedStocks(self, list_name, modelFunc, params, modelName=None, alarm='StoreDoor'):
        assert list_name in self.TrackLists
        trackList = self.TrackLists[list_name]
        stocksToTrack = trackList['stocksToTrack']
        stocksToTrackData = trackList['stocksToTrackData']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        stocksToTrackNextUpdateTime = trackList['stocksToTrackNextUpdateTime']
        stocksToTrackIndicators = trackList['stocksToTrackIndicators']

        for i in range(len(stocksToTrack)):
            symb =  stocksToTrack[i]
            df =  stocksToTrackData[i]
            if not isinstance(df, pd.DataFrame):
                continue
            if len(df)==0 or self.NpPdDatetimeToPyDatetime(df.index.values[-1]) < datetime.datetime.now() - datetime.timedelta(seconds=55):
                continue
            
            params_new = params.copy()
            for p in params:
                if params[p]==self.paramTypes.intradaylive:
                    params_new[p] = df
                if params[p]==self.paramTypes.intraday_prev:
                    assert False
                    params_new[p] = df_prev
                    
            if modelFunc(**params_new):
                print('='*60)
                print('{} hit: {}   @  {}'.format(modelName, symb, datetime.datetime.now()))
                print('='*60)
                self.SoundAlarm(alarm=alarm)
                
    
    
    
    ##remove inactive/delayed/no-data/very low volume/TODO: unhopeful stocks
    def cleanTrakedStocks(self, list_name, pv_30_th=100.0e3):
        assert list_name in self.TrackLists
        trackList = self.TrackLists[list_name]
        stocksToTrack = trackList['stocksToTrack']
        stocksToTrackData = trackList['stocksToTrackData']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']
        stocksToTrackNextUpdateTime = trackList['stocksToTrackNextUpdateTime']
        stocksToTrackIndicators = trackList['stocksToTrackIndicators']
        stocksRemovedFromTrackList = trackList['stocksRemovedFromTrackList']

        i = -1
        while True:
            i += 1
            if i>=len(stocksToTrack):
                break

            symb =  stocksToTrack[i]
            df =  stocksToTrackData[i]
            if not isinstance(df, pd.DataFrame):
                ##nodata
                self.RemoveStockFromTrackList(list_name, symb)
                continue
            if len(df)<1:
                ##nodata
                self.RemoveStockFromTrackList(list_name, symb)
                continue
            
            N = len(df)
            P = df['CLOSE'].values
            V = df['VOLUME'].values
            T = df.index.values
            t_n = T[-1]
            
            if t_n < np.datetime64(datetime.datetime.now()-datetime.timedelta(minutes=5)):
                ##delayed or inactive
                self.RemoveStockFromTrackList(list_name, symb)
                continue
            
            
            ##less than 100k in the last 30 minutes
            pv_30 = None
            dt = np.timedelta64(30, 'm')
            for j in range(N-2, 0, -1):
                if pv_30==None and (t_n-T[j])>=dt:
                    pv_30 = np.sum(P[j:]*V[j:])
                    break
            if pv_30!=None and pv_30<pv_30_th:
                self.RemoveStockFromTrackList(list_name, symb)
                continue
        
        return
    
        
    ##TODO: once a notable 1-2 min stock is noted monitor its real time proce for movements
    ##TODO: monitor notable energy companies every 10 mins
    ##TODO: monitor notable mining companies every 10 mins
    ## rep1

    def FillInRBCOrder(self, file_auth, password, account='TSFS:USD', symbol='AAPL', market='US', \
            order='buy', num_shares=1000, price=1000, AllorNone=True):
        from simplecrypt import encrypt, decrypt
        import os

        def read_encrypted(password, filename):
            with open(filename, 'rb') as input:
                ciphertext = input.read()
                plaintext = decrypt(password, ciphertext)
                return plaintext
        def decrypt_file(password, file_in, file_out):
            data_dec = read_encrypted(password, file_in)
            file_w = open(file_out, 'wb')
            file_w.write(data_dec)
            file_w.close()
            return True

        npz_name = 'tmp.npz'
        file_tmp = os.path.join(os.path.dirname(file_auth), npz_name)
        decrypt_file(password, file_auth, file_tmp)
        np_vars = np.load(file_tmp)
        user_ = str(np_vars['user_'])
        pass_ = str(np_vars['pass_'])
        #print(user_, '  ', pass_)
        os.remove(file_tmp)
        
        driver = webdriver.Chrome()
        driver.get(r"https://www1.royalbank.com/english/netaction/sgne.html")
        element = driver.find_element_by_id("USERID")
        element.clear()
        element.send_keys(user_)
        element = driver.find_element_by_id("PASSWORD")
        element.clear()
        element.send_keys(pass_)
        
        element = driver.find_element_by_xpath("//button[@title='Sign In']")
        element.click()
        WebDriverWait(driver, 10).until(EC.title_contains('Home - RBC Direct Investing'))
        
        #element = driver.find_element_by_css_selector('#accessible-megamenu-1488809184777-4.ul:nth-child(2).li:nth-child(2).a')
        #element.click()
        
        #element = driver.find_element_by_partial_link_text('Stocks &amp; ETFs')
        #element.click()
        
        #element = driver.find_element_by_xpath('//*[@id="accessible-megamenu-1488810487638-4"]/ul[1]/li[2]/a')
        #element.click()
        
        driver.get('https://www1.royalbank.com/cgi-bin/rbaccess/rbunxcgi?F22=4WN600S&7ASERVER=N601LD&LANGUAGE=ENGLISH&7ASCRIPT=/WebUI/OrderEntry/StockOrderEntryForm&r=BX5GdMiF_EOlHCg77uSB3w')
        
        #driver.find_element_by_xpath('//*[@id="selectedAction"]/option[3]').click()
        if order=='buy':
            driver.find_element_by_xpath("//select[@id='selectedAction']/option[@value='Buy']").click()
        elif order=='sell':
            driver.find_element_by_xpath("//select[@id='selectedAction']/option[@value='Sell']").click()
        
        element = driver.find_element_by_id("quantity")
        element.clear()
        element.send_keys(str(num_shares))

        element = driver.find_element_by_id("symbol")
        element.clear()
        element.send_keys(symbol)
        
        if market=='US':
            driver.find_element_by_xpath("//select[@id='selectedMarket']/option[@value='US']").click()
        if market=='CA':
            driver.find_element_by_xpath("//select[@id='selectedMarket']/option[@value='CA']").click()
        
        #time.sleep(2)
        driver.find_element_by_id("LimitPriceRadio").click()


        element = driver.find_element_by_id("LimitPriceValue")
        element.clear()
        element.send_keys(str(price))
        
        if AllorNone==True:
            #time.sleep(2)
            driver.find_element_by_id("AllOrNone").click()
        else:
            driver.find_element_by_id("AnyPartAccepted").click()

        driver.find_element_by_id("refreshQuote").click()

        return
        

    ##TODO: find latest target prices
    ##TODO: find latest Zack's ranking upgrades

            
    def TrimSortComapnyListsDF(self, trim=True, sort=True, hasSavedData=False, dir_path=None):
        if trim:
            self.df_NASDAQ = self.df_NASDAQ.loc[:,['Symbol','Name', 'MarketCap']]
            self.df_NYSE = self.df_NYSE.loc[:,['Symbol','Name', 'MarketCap']]
            self.df_NYSEMKT = self.df_NYSEMKT.loc[:,['Symbol','Name', 'MarketCap']]
        if sort:
            self.df_NASDAQ = self.df_NASDAQ.sort_values(by="MarketCap", ascending=False)
            self.df_NYSE = self.df_NYSE.sort_values(by="MarketCap", ascending=False)
            self.df_NYSEMKT = self.df_NYSEMKT.sort_values(by="MarketCap", ascending=False)
        if hasSavedData:
            assert dir_path!=None
            symbols = self.df_NASDAQ['Symbol']
            for i in range(len(symbols)):
                if not self.HasSavedData('NASDAQ:'+symbols[i], dir_path):
                    self.df_NASDAQ = self.df_NASDAQ[self.df_NASDAQ.Symbol != symbols[i]]
            symbols = self.df_NYSE['Symbol']
            for i in range(len(symbols)):
                if not self.HasSavedData('NYSE:'+symbols[i], dir_path):
                    self.df_NYSE = self.df_NYSE[self.df_NYSE.Symbol != symbols[i]]
            symbols = self.df_NYSEMKT['Symbol']
            for i in range(len(symbols)):
                if not self.HasSavedData('NYSEMKT:'+symbols[i], dir_path):
                    self.df_NYSEMKT = self.df_NYSEMKT[self.df_NYSEMKT.Symbol != symbols[i]]

        
    def HasSavedData(self, symb, dir_path):
        file_name = symb.replace(':', '_')
        file_path = os.path.join(dir_path, file_name)
        if os.path.exists(file_path):
            return True
        else:
            return False
            
    def NpPdDatetimeToPyDatetime(self, date_time):
        if isinstance(date_time, np.datetime64):
            return date_time.astype('M8[ms]').astype('O')
        elif isinstance(date_time, pd.tslib.Timestamp):
            return date_time.to_pydatetime()
        elif isinstance(date_time, datetime.datetime):
            return date_time
    
    def NpPdtimedeltaToPytimedelta(self, dt):
        if isinstance(dt, np.timedelta64):
            return datetime.timedelta(seconds=dt.astype('timedelta64[ms]').item().total_seconds())
        elif isinstance(dt, pd.Timedelta):
            return dt.to_pytimedelta()
        elif isinstance(dt, datetime.timedelta):
            return dt

    def NpPdtimedeltaArrayToMins(self, dT):
        N = len(dT)
        dT_min = np.zeros(N)
        for i in range(N):
            dT_min[i] = int(self.NpPdtimedeltaToPytimedelta(dT[i]).seconds/60)
        return dT_min
    

            
    def UpdateIntradayStockData(self, dir_path, market='NASDAQ', interval_sec=60, period_days=30):
        assert market in ['NASDAQ', 'NYSE', 'NYSEMKT']
        
        company_list = None
        if market=='NASDAQ':
            company_list = self.df_NASDAQ
        elif market=='NYSE':
            company_list = self.df_NYSE
        elif market=='NYSEMKT':
            company_list = self.df_NYSEMKT
            
        assert 'Symbol' in company_list.keys()
        
        company_symbols = company_list['Symbol']
        n_comp = len(company_symbols)
        
        time_last = time.perf_counter()
        dowload_sec = 29
        stop_sec = 31
        
        for i in range(len(company_symbols)):
            if time.perf_counter()-time_last>dowload_sec:
                time.sleep(stop_sec)
                time_last = time.perf_counter()
        
            s = company_symbols[i]
            symb = s.strip()
            mkt_symb = market+':'+symb
            print('{}:{}/{}'.format(mkt_symb, i, n_comp), end=' ')
            #print(symb, end=' ')
            file_name = mkt_symb.replace(':', '_')

            file_path  = os.path.join(dir_path, file_name)
            
            if not os.path.exists(file_path):
                endDate = datetime.date.today()
                
                try:
                    df = self.getIntradayOrHistoricPrices(mkt_symb, str(interval_sec), str(period_days)+'d')
                    if len(df)>0:
                        df = df.set_index('DATE')
                        
                        df.to_hdf(file_path, symb, format='table', mode='w')
                except Exception as e:
                    print(symb, ' error: ', str(e))
                except:
                    print(symb, ' error: ', sys.exc_info()[0])
            else:
                store = pd.HDFStore(file_path)
                keys = store.keys()
                if r'/'+symb not in keys:
                    os.remove(file_path)
                    print('\n {} deleted.'.format(file_path))
                    continue
                nrows = store.get_storer(symb).nrows
                
                start = None
                if nrows>10:
                    start = nrows-10
                #df = pd.read_hdf(file_path, symb, start=start)
                df = store.select(symb, start=start)
                store.close()
                
                df.sort_index(ascending=True,  inplace=True)
                period = period_days
                if len(df)>0:
                    date_last = self.NpPdDatetimeToPyDatetime(df.index[0])
                    date_today = datetime.datetime.today()
                    period = (date_today - date_last).days + 1
                
                
                period = min(period, period_days)
                period = str(period)+'d'

                try:
                    df = self.getIntradayOrHistoricPrices(mkt_symb, str(interval_sec), period)
                    df = df[df.DATE > date_last]
                    df = df.set_index('DATE')
                    
                    df.to_hdf(file_path, symb, format='table', mode='a', append=True)
                except Exception as e:
                    print(symb, ' error: ', repr(e))
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname)
                    traceback.print_exc()
                except:
                    print(symb, ' error: ', sys.exc_info()[0])

    

    def CorrectMarketDataErrors(self, df):
        ## Drop duplicates
        df = df.reset_index().drop_duplicates(subset='DATE', keep='first').set_index('DATE')
        
        ##remove errors
        T = df.index.values
        N = len(T)
        ind_error = []
        
        for i in range(1,N-1):
            if not (T[i]>T[i-1] and T[i]<T[i+1]):
                ind_error.append(i)
                
        df = df.drop(df.index[ind_error])
        return df
    

    def UpdateIntradayStockData_All(self):
        self.UpdateIntradayStockData(self.companyDataIntraday_path, market='NASDAQ')
        self.UpdateIntradayStockData(self.companyDataIntraday_path, market='NYSE')
        self.UpdateIntradayStockData(self.companyDataIntraday_path, market='NYSEMKT')

    def UpdateIntraday30minStockData_All(self):
        self.UpdateIntradayStockData(self.companyDataIntraday30min_path, market='NASDAQ', interval_sec=1800, period_days=500)
        self.UpdateIntradayStockData(self.companyDataIntraday30min_path, market='NYSE', interval_sec=1800, period_days=500)
        self.UpdateIntradayStockData(self.companyDataIntraday30min_path, market='NYSEMKT', interval_sec=1800, period_days=500)


    def plotSymbolIntraday(self, symb, duration=[5, 'd'], curves=['CLOSE', 'VOLUME'], window=None, figsize=(10, 6),
                    extplt=False, extpltData=[10, 1]):
        """ extpltData = [periode (days), number of extrapolations]
        """
        file_name = symb.replace(':', '_')
        file_path  = os.path.join(self.companyDataIntraday_path, file_name)
        if os.path.exists(file_path):
            assert symb.find(':')>=0
            symb_bare = symb[symb.find(':')+1:]
            df = pd.read_hdf(file_path, symb_bare)
            
            df = self.CorrectMarketDataErrors(df)
            
            if window!=None:
                df = df.rolling(window=window).mean()
            
            len_mins = len(df)
            if duration[1]=='Y':
                len_mins = int(duration[0]*365*8*60)
            if duration[1]=='d':
                len_mins = int(duration[0]*8*60)
                
            if len_mins<len(df):
                df = df[-len_mins:]

            n_p = len(curves)
            f, axarr = plt.subplots(n_p, sharex=True, figsize=figsize)

            p_i = 0
            if 'CLOSE' in curves:
                axarr[p_i].plot(df['CLOSE'].values)
                N_df = len(df)
                if extplt:
                    mins_extplt, n_extplt = extpltData
                    for i in range(n_extplt):
                        y = df['CLOSE']
                        y = y[np.logical_not(np.isnan(y))]
                        n_y = (i+1)*mins_extplt
                        if n_y<len(y):
                            y_ = y[-n_y:]
                            x_ = np.arange(len(y_))
                            y_L = polynomial.Polynomial(polynomial.polyfit(x_, y_, 1))(x_) 
                            axarr[p_i].plot(np.arange(N_df)[-len(y_L):], y_L)
                p_i += 1
            if 'VOLUME' in curves:
                axarr[p_i].plot(df['VOLUME'].values/1.e6)
                p_i += 1
            plt.show()
            return df
        else:
            print('Market data for {} not available.'.format(symb))
            return None
            
    def GetIntradayMinutData(self, symb, day):
        dir_path = self.companyDataIntraday_path
        file_name = symb.replace(':', '_')
        file_path  = os.path.join(dir_path, file_name)
    
        assert symb.find(':')>=0
        symb_bare = symb[symb.find(':')+1:]

        if os.path.exists(file_path):
            store = pd.HDFStore(file_path)
            keys = store.keys()
            if r'/'+symb_bare not in keys:
                return None
            nrows = store.get_storer(symb_bare).nrows
            
            df = store.select(symb_bare, [pd.Term('index', '>=', day)])
            store.close()
            
            df = self.CorrectMarketDataErrors(df)
            
            day = datetime.datetime.combine(day, datetime.time(hour=9, minute=30))
            day_open  = day.replace(hour=9, minute=30, second=0, microsecond=0)
            day_close = day.replace(hour=16, minute=0, second=0, microsecond=0)

            df = df[np.logical_and(df.index>=day_open, df.index<=day_close)]
            
            return df



    def FindNoticableIntradayHistoricalSymbols(self, period_days=30, perc_change=10, mins_change=None, GL=None):
        """ GL = 'G': gainers,   'L': losers    None: both
        changing perc_change percent in mins_change minutes
        """
        dir_path = self.companyDataIntraday_path

        ## find thos companies whose prices has changed by perc_change in the last period_days 
        markets = ['NASDAQ', 'NYSE', 'NYSEMKT']
        
        noticable_dic = {}
        day_start = datetime.date.today() - datetime.timedelta(days=period_days)
        for market in markets:
            company_list = None
            if market=='NASDAQ':
                company_list = self.df_NASDAQ
            elif market=='NYSE':
                company_list = self.df_NYSE
            elif market=='NYSEMKT':
                company_list = self.df_NYSEMKT
                
            assert 'Symbol' in company_list.keys()
            
            company_symbols = company_list['Symbol']
            n_comp = len(company_symbols)
            
            for ind_comp in range(len(company_symbols)):
                s = company_symbols[ind_comp]
                symb = s.strip()
                mkt_symb = market+':'+symb
                file_name = mkt_symb.replace(':', '_')
                file_path  = os.path.join(dir_path, file_name)
                
                if os.path.exists(file_path):
                    store = pd.HDFStore(file_path)
                    keys = store.keys()
                    if r'/'+symb not in keys:
                        continue
                    nrows = store.get_storer(symb).nrows
                    
                    df = store.select(symb, [pd.Term('index', '>', day_start)])
                    store.close()
                    
                    df = self.CorrectMarketDataErrors(df)
                    for del_day in range(period_days):
                        day_j = datetime.datetime.combine(day_start + datetime.timedelta(days=del_day), datetime.time(hour=9, minute=30))
                        day_open  = day_j.replace(hour=9, minute=30, second=0, microsecond=0)
                        day_close = day_j.replace(hour=16, minute=0, second=0, microsecond=0)

                        df_j = df[np.logical_and(df.index>=day_open, df.index<=day_close)]
                        
                        if len(df_j)>100:
                            P = df_j['CLOSE'].values
                            T = df_j.index.values
                            #print(len(P), len(T))
                            assert len(P)==len(T)
                            
                            if (P.max()-P.min())/P.mean()*100 > perc_change:
                                if mins_change==None:
                                    ind_max = np.argmax(P)
                                    ind_min = np.argmin(P)
                                    if GL=='GL' or GL==None:
                                        if mkt_symb in noticable_dic:
                                            noticable_dic[mkt_symb].append(day_open.date())
                                        else:
                                            noticable_dic[mkt_symb] = [day_open.date()]
                                    elif GL=='G' and ind_max>ind_min:
                                        if mkt_symb in noticable_dic:
                                            noticable_dic[mkt_symb].append(day_open.date())
                                        else:
                                            noticable_dic[mkt_symb] = [day_open.date()]
                                    elif GL=='L' and ind_max<ind_min:
                                        if mkt_symb in noticable_dic:
                                            noticable_dic[mkt_symb].append(day_open.date())
                                        else:
                                            noticable_dic[mkt_symb] = [day_open.date()]
                                else:
                                    T0 = T-T[0]
                                    T0_mins = np.zeros(len(T0))
                                    for i in range(len(T0)):
                                        T0_mins[i] = self.NpPdtimedeltaToPytimedelta(T0[i]).seconds//60
                                    passed = False
                                    gained = False
                                    for i in range(len(T)-mins_change-1):
                                        if passed:
                                            break
                                        for j in range(1, mins_change+1):
                                            if T0_mins[i+j]-T0_mins[i]>mins_change:
                                                break
                                            if abs(P[i+j]-P[i])/P[i]*100>perc_change:
                                                passed = True
                                                if P[i+j]>P[i]:
                                                    gained = True
                                                print(mkt_symb, T[i], P[i], P[i+j])
                                                break
                                    if passed:
                                        if GL=='GL' or GL==None:
                                            if mkt_symb in noticable_dic:
                                                noticable_dic[mkt_symb].append(day_open.date())
                                            else:
                                                noticable_dic[mkt_symb] = [day_open.date()]
                                        elif GL=='G' and gained:
                                            if mkt_symb in noticable_dic:
                                                noticable_dic[mkt_symb].append(day_open.date())
                                            else:
                                                noticable_dic[mkt_symb] = [day_open.date()]
                                        elif GL=='L' and not gained:
                                            if mkt_symb in noticable_dic:
                                                noticable_dic[mkt_symb].append(day_open.date())
                                            else:
                                                noticable_dic[mkt_symb] = [day_open.date()]
                                                                     
                                
        return noticable_dic                                


    ## Models
    """ Model 1: volume breaks last day, for first time, it is the highest intraday volume, 
        and price increased in the last minute more than the open jump
        name: ModelVolBreakPriceup
    """
    def ModelVolBreakPriceup(self, df, df_prev, period_days=30, vol_prev_ratio=1.0, vol_intra_ratio=2.0, 
                price_open_jump_max=10, price_jump_1min=2, price_min=1.0, price_max=10.0):
        df_j = df
        
        V = df_j['VOLUME'].values
        P = df_j['CLOSE'].values
        T = df_j.index.values
        
        if P.min()<price_min or P.max()>price_max:
            return False

        dT = T[1:] - T[:-1]
        dT_mins = self.NpPdtimedeltaArrayToMins(dT)
        
        df_j_prev = df_prev
        _V = df_j_prev['VOLUME'].values
        _P = df_j_prev['CLOSE'].values
        _T = df_j_prev.index.values

        if len(_V)<100:
            return False

        _P_close = _P[-1]
        _V_max = _V.max()
        
        P_open_jump = (P[0]-_P_close)/_P_close*100
        if P_open_jump>price_open_jump_max:
            return False

        N = len(T)
        i = N-1
        vol_break = False                            
        if V[i]>_V_max*vol_prev_ratio:
            vol_break = True

        price_up = False
        if (P[i]-P[i-1])/P[i-1]*100>price_jump_1min:
            price_up = True

        vol_intra_break = False
        if (V[i]-V[:i].max())/V[:i].max()*100>vol_intra_ratio:
            vol_intra_break = True
        
        if vol_break and price_up and vol_intra_break:
            return True

        return False


    def ModelintVdPSharpFall(self, df):
        """ The proce dropes sharply
        """
        if len(df)<2:
            return False
        P = df['CLOSE'].values
        V = df['VOLUME'].values
        T = df.index.values
        
        N = len(T)

        dP = P[1:]-P[:-1]
        dV = V[1:]-V[:-1]
        dT = T[1:]-T[:-1]
        
        N_avg = 5
        P_avg = np.zeros(N)
        for i in range(N):
            P_avg[i] = np.sum(P[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))

        N_avg = 10
        V_avg = np.zeros(N)
        for i in range(N):
            V_avg[i] = np.sum(V[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))

        dPV = P[1:]*V[1:]-P[:-1]*V[:-1]
        VdP = V[1:]*dP
        intVdP = np.zeros(len(VdP))
        intVdP[0] = VdP[0]
        for i in range(1, len(VdP)):
            intVdP[i] = intVdP[i-1]+VdP[i]
        intVdP_rel = intVdP/np.abs(intVdP).max()

        dintVdP_rel = intVdP_rel[1:] - intVdP_rel[:-1]

        N_avg = 20
        intVdP_avg = np.zeros(N-1)
        for i in range(N-1):
            intVdP_avg[i] = np.sum(intVdP[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))
        intVdP_rel_avg = intVdP_avg/np.abs(intVdP).max()
        
        i_up = -10
        N = len(T)
        hit_k , miss_k = 0, 0
        hit_g , miss_g = 0, 0
        i = N-3
        if np.abs(P[:i+2]).max() - np.abs(P[:i+2]).min()>np.abs(P[:i+2]).min()/10:
            if dintVdP_rel[i]<-np.abs(dintVdP_rel[:i+1]).max()/5 and \
                (P[i+2]-P_avg[i+2])/P_avg[i+2]<-0.05:
                if np.argmax(intVdP_rel[:i+1])<i-5:
                    return True
        return False
         

    def ModelintVdPStopThenIncrease(self, df, p_change_min=0.05, p_change_max=0.15, ind_max_prev=3, confirm_incease=False):
        """ The price velocity decreases momentarily and it continues to increase
        """
        if len(df)<2:
            return False
        P = df['CLOSE'].values
        V = df['VOLUME'].values
        T = df.index.values
        
        N = len(T)
        if N<2:
            return False

        dP = P[1:]-P[:-1]
        dV = V[1:]-V[:-1]
        dT = T[1:]-T[:-1]
        
        N_avg = 5
        P_avg = np.zeros(N)
        for i in range(N):
            P_avg[i] = np.sum(P[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))

        N_avg = 10
        V_avg = np.zeros(N)
        for i in range(N):
            V_avg[i] = np.sum(V[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))

        dPV = P[1:]*V[1:]-P[:-1]*V[:-1]
        VdP = V[1:]*dP
        intVdP = np.zeros(N)
        intVdP[0] = VdP[0]
        for i in range(1, N):
            intVdP[i] = intVdP[i-1]+VdP[i-1]
        intVdP_rel = intVdP/np.abs(intVdP).max()

        dintVdP_rel = intVdP_rel[1:] - intVdP_rel[:-1]

        N_avg = 20
        intVdP_avg = np.zeros(N-1)
        for i in range(N-1):
            intVdP_avg[i] = np.sum(intVdP[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))
        intVdP_rel_avg = intVdP_avg/np.abs(intVdP).max()
        
        i = N-1
        if (P[max(i-10,0):].max()-P[max(i-10,0):].min())/P[max(i-10,0):].min()>p_change_min and\
            (P[max(i-10,0):].max()-P[max(i-10,0):].min())/P[max(i-10,0):].min()<p_change_max and np.argmax(P)>i-ind_max_prev:
            ind_max = np.argmax(intVdP_rel)
            if max(i-10,2)<ind_max<=i:
                if confirm_incease:
                    if np.argmin(intVdP_rel[ind_max:])<(i+1-ind_max) and np.argmax(intVdP_rel)==i:
                        return True
                    elif np.argmax(intVdP_rel)==i and intVdP_rel[i-1]>0 and dintVdP_rel[i]>0.2*dintVdP_rel[i-1]:
                        return True
                elif np.argmax(intVdP_rel)==i-1:
                    return True
                
        return False
    

    def ModelSharpRise(self, df, p_change_min=0.05, price_vol_min=50.0e3, price_min=1.0, price_max=10.0):
        """ The price velocity decreases momentarily and it continues to increase
        """
        if len(df)<2:
            return False
        P = df['CLOSE'].values
        V = df['VOLUME'].values
        T = df.index.values
        
        N = len(T)
        if N<2:
            return False

        dP = P[1:]-P[:-1]
        dV = V[1:]-V[:-1]
        dT = T[1:]-T[:-1]
        
        N_avg = 5
        P_avg = np.zeros(N)
        for i in range(N):
            P_avg[i] = np.sum(P[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))

        N_avg = 10
        V_avg = np.zeros(N)
        for i in range(N):
            V_avg[i] = np.sum(V[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))

        dPV = P[1:]*V[1:]-P[:-1]*V[:-1]
        VdP = V[1:]*dP
        intVdP = np.zeros(N)
        intVdP[0] = VdP[0]
        for i in range(1, N):
            intVdP[i] = intVdP[i-1]+VdP[i-1]
        intVdP_rel = intVdP/np.abs(intVdP).max()

        dintVdP_rel = intVdP_rel[1:] - intVdP_rel[:-1]

        N_avg = 20
        intVdP_avg = np.zeros(N-1)
        for i in range(N-1):
            intVdP_avg[i] = np.sum(intVdP[max(i+1-N_avg, 0):i+1])/(i+1-max(i+1-N_avg, 0))
        intVdP_rel_avg = intVdP_avg/np.abs(intVdP).max()
        
        i_rise = -30
        i = N-1
        for i in range(2, N):
            ii = i-1
            if np.abs(P[:i+1]).max() - np.abs(P[:i+1]).min()>np.abs(P[:i+1]).min()*p_change_min \
                    and P[i]*V[i]>price_vol_min and price_min<P[i]<price_max and V[i]>0.9*V[i-1]:
                X = np.abs(intVdP_rel[:max(ii-2,2)]*(intVdP_rel[:max(ii-2,2)]>0.0)).max()
                if (i-i_rise)>25 and intVdP_rel[ii]>2.0*X:# and intVdP_rel[ii]<10.0*X:
                    i_rise = i
                    if i_rise==N-1:
                        return True
        return False


    def ModelPriceOnly_SharpRise(self, df, price_change_percent=10.0, time_change=1):
        if len(df)<2:
            return False
        P = df['CLOSE'].values
        V = df['VOLUME'].values
        T = df.index.values
        
        N = len(T)
        if N<2:
            return False
            
        if P[-1]-P[-2]>price_change_percent*P[-2]/100.0:
            return True
            
        return False
     



    def TestModels(self, modelFunc, params, symbols_dates=None, day_start=None, \
            day_end=None, price_change_perc=[3, 10]):
        """ price_change_perc[3, 10] : 3 percent in 10 minutes
        """
        dir_path = self.companyDataIntraday_path
        if symbols_dates==None:
            assert day_start!=None and day_end!=None
            symbols_dates = {}
            
            dt_list = pd.bdate_range(start=day_start, end=day_end)
            dates_list = dt_list.date.tolist() 
            
            markets = ['NASDAQ', 'NYSE', 'NYSEMKT']
            for market in markets:
                company_list = None
                if market=='NASDAQ':
                    company_list = self.df_NASDAQ
                elif market=='NYSE':
                    company_list = self.df_NYSE
                elif market=='NYSEMKT':
                    company_list = self.df_NYSEMKT
                    
                assert 'Symbol' in company_list.keys()
                
                company_symbols = company_list['Symbol']
                n_comp = len(company_symbols)
                
                symbols_dates.update(dict(zip(company_symbols, dates_list)))
                
        n_hit, n_miss = 0, 0
        for mkt_symb in symbols_dates:
            assert mkt_symb.find(':')>=0
            symb = mkt_symb[mkt_symb.find(':')+1:]

            file_name = mkt_symb.replace(':', '_')
            file_path  = os.path.join(dir_path, file_name)
            
            dates = sorted(symbols_dates[mkt_symb])
            
            if os.path.exists(file_path):
                store = pd.HDFStore(file_path)
                keys = store.keys()
                if r'/'+symb not in keys:
                    continue
                nrows = store.get_storer(symb).nrows
                
                df = store.select(symb, [pd.Term('index', '>=', dates[0]-BDay(1))])
                store.close()
                
                df = self.CorrectMarketDataErrors(df)
                for day in dates:
                    day_open = datetime.datetime.combine(day, datetime.time(hour=9, minute=30))
                    day_close = day_open.replace(hour=16, minute=0, second=0, microsecond=0)

                    df_j = df[np.logical_and(df.index>=day_open, df.index<=day_close)]
                    
                    if len(df_j)<100:
                        continue
                        
                    day_prev = day-BDay(1)
                    day_prev_open = datetime.datetime.combine(day_prev, datetime.time(hour=9, minute=30))
                    day_prev_close = day_prev_open.replace(hour=16, minute=0, second=0, microsecond=0)
                        
                    df_j_prev = df[np.logical_and(df.index>=day_prev_open, df.index<=day_prev_close)]
                    
                    if self.paramTypes.intraday_prev in list(params.values()):
                        if len(df_j_prev)<100:
                            continue
                    
                    for i in range(len(df_j)):
                    
                        params_new = params.copy()
                        for p in params:
                            if params[p]==self.paramTypes.intradaylive:
                                params_new[p] = df_j[:i+1]
                            if params[p]==self.paramTypes.intraday_prev:
                                params_new[p] = df_j_prev
                                
                        if modelFunc(**params_new):
                            P = df['CLOSE']
                            dp_p, t_min = price_change_perc
                            if (P[i:max(i+t_min, len(df_j))].max()-P[i])/P[i]*100>=dp_p:
                                n_hit += 1
                            else:
                                n_miss += 1
        print('n_hit: {}   n_miss: {}'.format(n_hit, n_miss))
        if n_miss>0:
            print('hit/miss ratio: {}'.format(n_hit/n_miss))
        


    def FindPremarketMostActives(self, price_min=0.9, price_max=2.5, n_trans_min=20):
        dir_path = self.companyDataIntraday_path
        active_premarkets = []
        
        t_start = time.perf_counter()
        
        markets = ['NASDAQ', 'NYSE', 'NYSEMKT']
        for market in markets:
            company_list = None
            if market=='NASDAQ':
                company_list = self.df_NASDAQ
            elif market=='NYSE':
                company_list = self.df_NYSE
            elif market=='NYSEMKT':
                company_list = self.df_NYSEMKT
                
            assert 'Symbol' in company_list.keys()
            
            df = company_list
            df = df[df['LastSale']!='n/a']
            df['LastSale'] = df['LastSale'].astype(float).fillna(0.0)
            df = df[np.logical_and(df['LastSale']<price_max, df['LastSale']>price_min)]

            company_symbols = df['Symbol']
            n_comp = len(company_symbols)
            df = None

            time_last = time.perf_counter()
            dowload_sec = 29
            stop_sec = 31

            for symb_bare in company_symbols:
                if time.perf_counter()-time_last>dowload_sec:
                    time.sleep(stop_sec)
                    time_last = time.perf_counter()
    
                symb = market+':'+symb_bare

                period = '2d'
                intval = '60'

                df = None
                new_data_available = False
                try:
                    df = self.getIntradayOrHistoricPrices(symb, intval, period, extendex_hours=True)
                    new_data_available = True
                except:
                    pass
                if new_data_available:
                    assert isinstance(df, pd.DataFrame)
                    assert 'DATE' in df
                    
                    df = df.set_index('DATE')
                    #print(df, len(df), end=' ')
                    
                    t_now = datetime.datetime.now()
                    t_open = t_now.replace(hour=9, minute=30, second=0, microsecond=0)
                    t_earlymorning = t_now.replace(hour=6, minute=0, second=0, microsecond=0)
                    df = df[np.logical_and(df.index<t_open, df.index>t_earlymorning)]
                
                    #print('{}-({}) '.format(symb, len(df)), end=' ')

                    if len(df)>n_trans_min:
                        active_premarkets.append([symb, df])
                        V = df['VOLUME'].values
                        P = df['CLOSE'].values
                        print('{} : N: {}    V: {}k    P: {}k'.format(symb, len(df), int(V.sum()/1.0e3), int((V*P).sum()/1.0e3)))
                            
        t_end = time.perf_counter()
        print('total time: ', (t_end-t_start)/60)

        return active_premarkets                            


    def FindAftermarketMostActives(self, price_min=0.9, price_max=2.5, n_trans_min=20):
        dir_path = self.companyDataIntraday_path
        active_aftermarkets = []
        
        t_start = time.perf_counter()
        
        markets = ['NASDAQ', 'NYSE', 'NYSEMKT']
        for market in markets:
            company_list = None
            if market=='NASDAQ':
                company_list = self.df_NASDAQ
            elif market=='NYSE':
                company_list = self.df_NYSE
            elif market=='NYSEMKT':
                company_list = self.df_NYSEMKT
                
            assert 'Symbol' in company_list.keys()
            
            df = company_list
            df = df[df['LastSale']!='n/a']
            df['LastSale'] = df['LastSale'].astype(float).fillna(0.0)
            df = df[np.logical_and(df['LastSale']<price_max, df['LastSale']>price_min)]

            company_symbols = df['Symbol']
            n_comp = len(company_symbols)
            df = None

            time_last = time.perf_counter()
            dowload_sec = 29
            stop_sec = 31

            for symb_bare in company_symbols:
                if time.perf_counter()-time_last>dowload_sec:
                    time.sleep(stop_sec)
                    time_last = time.perf_counter()
    
                symb = market+':'+symb_bare

                period = '2d'
                intval = '60'

                df = None
                new_data_available = False
                try:
                    df = self.getIntradayOrHistoricPrices(symb, intval, period, extendex_hours=True)
                    new_data_available = True
                except:
                    pass
                if new_data_available:
                    assert isinstance(df, pd.DataFrame)
                    assert 'DATE' in df
                    
                    df = df.set_index('DATE')
                    #print(df, len(df), end=' ')
                    
                    t_now = datetime.datetime.now()
                    t_close = t_now.replace(hour=16, minute=00, second=0, microsecond=0)
                    t_night = t_now.replace(hour=20, minute=0, second=0, microsecond=0)
                    df = df[np.logical_and(df.index>t_close, df.index<t_night)]
                
                    #print('{}-({}) '.format(symb, len(df)), end=' ')

                    if len(df)>n_trans_min:
                        active_aftermarkets.append([symb, df])
                        V = df['VOLUME'].values
                        P = df['CLOSE'].values
                        print('{} : N: {}    V: {}k    P: {}k'.format(symb, len(df), int(V.sum()/1.0e3), int((V*P).sum()/1.0e3)))
                            
        t_end = time.perf_counter()
        print('total time: ', (t_end-t_start)/60)

        return active_aftermarkets                            


    def SetLatestRecordedPrices(self):
        
        return

        
##------------------------------------------------------------------        
        
class StockReports:
    def __init__(self, form_type='8-K', form_items=['2.02'], count=100):
        """
            form_type = '10-Q', form_items = None
            form_type = '8-K',  form_items = ['2.02']#, '9.01']
        """
        self.SECLink = 'https://www.sec.gov'
        
        self.form_type = form_type
        self.form_items = form_items
        self.count = count
        
        self.vbose = False

        self.mp3alarm_path = os.path.join('other', 'Store_Door_Chime.mp3')

        self.yeartitles = ['2015', '2016']
        self.quartertitles = ['1Q15', '2Q15', '3Q15', '4Q15', '1Q16', '2Q16', '3Q16', '4Q16']
        
        self.Y_orders = dict(zip(self.yeartitles, range(len(self.yeartitles))))
        self.Q_orders = dict(zip(self.quartertitles, range(len(self.quartertitles))))
        
        ## '3MY': three months yearly,  'YY': no three/nine month label
        class ColLabels(Enum):
            M3Y = 0
            M6Y = 1
            M9Y = 2
            M12Y = 3
            YY = 4 
            W39Y = 5
            W13Y = 6
            Q1Y = 7
            Q2Y = 8
            Q3Y = 9
            Q4Y = 10

        self.ColLabels = ColLabels

        ## description labels
        class descLabels(Enum):
            totalAssets = 0
            netProfit = 1
            totalLiabilities = 2
            netProfitPerShareBasic = 3 
            netProfitPerShareDiluted = 4 
            EBITDA = 5


        self.descLabels = descLabels
        
        self.logfile = os.path.join('Finance', 'log.txt')
        
        # _pm : initial value positive, change negative(minus)
        self.weightingFactors_pp = {'asset':1.0,  'liability':-1.0, 'profit':4, 'profit per share':4}
        self.weightingFactors_pm = {'asset':1.0,  'liability':-1.0, 'profit':4, 'profit per share':4}
        self.weightingFactors_mp = {'asset':None, 'liability':None, 'profit':-3, 'profit per share':-3}
        self.weightingFactors_mm = {'asset':None, 'liability':None, 'profit':+3, 'profit per share':+3}
        
        return
        

    def logToFile(self, txtLine):
        with open(self.logfile, 'a') as file_:
            file_.write('\n' + txtLine + '\n')                
            file_.write('-'*50 + '\n')                
    
    
    def GetSecRecetlyReportedFeed(self, page_no=0):
        form_type = self.form_type
        count = self.count
        start = page_no*count

        edgar_recent_rss_url = 'https://www.sec.gov/cgi-bin/browse-edgar?action=getcurrent&CIK=&type={}&company=&dateb=&owner=include&start={}&count={}&output=atom'.format(form_type, start, count) 
        feed = feedparser.parse( edgar_recent_rss_url )
        
        self.feed = feed
        return feed


    def PrintFeedInfo(self, feed):
        print('feed[ "bozo" ]: ', feed[ "bozo" ])
        print('feed[ "url" ]: ', feed[ "url" ])
        print('feed[ "version" ]: ', feed[ "version" ])
        print('feed[ "channel" ]: ', feed[ "channel" ].keys() )
        for i in range(len(feed[ "items" ])):
            print(i, ' : ', feed["items"][i]['title'], feed["items"][i]['summary_detail']['value'], '\n')
            
        if len(feed[ "items" ])>0:
            print('feed[ "items" ]: ', feed[ "items" ][0].keys(), '\n ', feed[ "items" ][0])


    def FindFeeditemsHavingFormitem(self):
        if self.form_items==None or len(self.form_items)==0:
            return list(range(self.count))
        feed = self.feed
        form_items = self.form_items
        forms_with_items = []
        for i in range(self.count):
            has_form_item = False
            if form_items==None or len(form_items)==0:
                has_form_item = True
            for j in range(len(form_items)):
                if str(feed["items"][i]).find('Item '+form_items[j])>=0:
                    has_form_item = True
                    forms_with_items.append(i)
                    break

        if self.vbose:
            print('Forms with requested items: ', forms_with_items)

        return forms_with_items
        


    #TODO: process table to structured data
    """
    Find blocks (continuous set of rows) with the same number of columns
    find if the block has numeric information.. like numbers in the same  column
    what comes before is title
    insid the block find the subtables
    extract subtable data
    drop $ % )..
    drop/check single closing parenthesis
    """

    def HtmlTable_ExtractData(self, html_table):
        ##TODO
        rows = table.findAll(lambda tag: tag.name=='tr')
        n_col_max = 0
        n_cols_last = None
        block_r_start = 0
        
        for i in range(len(rows)-1,-1,-1):
            row = rows[i]
            cols = row.findAll('td')
            n_col = 0
            for col in cols:
                if col.has_attr('colspan'):
                    n_col += int(col['colspan'])
                else:
                    n_col += 1
            if n_col>n_col_max:
                n_col_max = n_col
            if n_cols_last==None:
                n_cols_last = n_col
            if n_col!=n_cols_last:
                block_r_start = i+1
        if block_r_start<len(rows)/2:
            ##get title info
            description = ""
            if block_r_start>0:
                for i in range(block_r_start):
                    row = rows[i]
                    cols = row.findAll('td')
                    n_col = 0
                    for col in cols:
                        if col.renderContents().find('<div')>=0:
                            print('Has div..')
                        else:
                            description += col.renderContents()
            print('description: ', description)
            ##
        else:
            print('Unrecognized table structure.')
            return None



    def isNumeric(self, expr):
        for c in '.,()$%':
            expr = expr.replace(c,'')
        return expr.isdigit()

    def isDescription(self, string):
        if not re.search('[e-zE-Z]', string):
            return False
        else:
            return True
            
    def numericStrToFloat(self, _expr_):
        expr = _expr_
        for c in ')$%':
            expr = expr.replace(c,'')
        expr = expr.replace('(','-')
        ind_r_comma = expr[::-1].find(',')
        if ind_r_comma==3:
            expr = expr.replace(',','')
        elif ind_r_comma>=0 and '.' not in expr:
            n_comma = 0
            for c in expr:
                if c==',':
                    n_comma += 1
            if n_comma<=1:
                expr = expr.replace(',', '.')
            else:
                print('Number format {} not recognized.'.format(_expr_))
                self.logToFile('Number format {} not recognized.'.format(_expr_))
                return None
        elif ind_r_comma>=0 or '.' in expr:
            expr = expr.replace(',','')

        try:
            return float(expr)
        except:
            print('Error: float({}) not valid.'.format(expr))
            return None
            
        

    def isColTitle(self, string):
        for title in self.yeartitles:
            if title in string:
                return True
        for title in self.quartertitles:
            if title in string:
                return True
        return False

    """
    def GetColTitleLabel(self, title):
        ##TODO: treat quarter labels as well
        title = title.lower()
        for t in self.yeartitles:
            if t in label:
                if 'three month' in title or '3 month' in title or 'quarter ended' in title:
                    return self.ColLabels.M3Y, self.Y_orders[t]
                if 'nine month' in title or '9 month' in title:
                    return self.ColLabels.M9Y, self.Y_orders[t]
                if 'thirteen weeks' in title or '13 weeks' in title or\
                   'fourteen weeks' in title or '14 weeks' in title or\
                   'twelve weeks' in title or '12 weeks' in title:
                    return self.ColLabels.W13Y, self.Y_orders[t]
                if 'thirty-nine weeks' in title or 'thirty nine weeks' in title or '39 weeks' in title or\
                   'thirty-eight weeks' in title or 'thirty eight weeks' in title or '38 weeks' in title or\
                   'fourty weeks' in title or '40 weeks' in title:
                    return self.ColLabels.W39Y, self.Y_orders[t]
                if 'quarter' not in title:
                    return self.ColLabels.YY, self.Y_orders[t]

        return None
    """
        
    def HtmlTable_GetTermAndPeriod(self, table, term='net profit'):
        if str(table).lower().find(term)<0:
            return None
        rows = table.findAll(lambda tag: tag.name=='tr')
        n_col_max = 0
        n_cols_last = None
        block_r_start = 0
        
        financialData = []
        for i in range(len(rows)-1,-1,-1):
            row = rows[i]
            if str(row).lower().find(term)<0:
                continue
            cols = row.findAll('td')
            n_col = 0
            data_cols = []
            for col in cols:
                #print(n_col, end=':')
                c_ind_st = n_col
                c_len = 1
                if col.has_attr('colspan'):
                    c_len = int(col['colspan'])
                n_col += c_len
                col_text = col.text.strip()
                #print(col_text, end='  |')
                if self.isNumeric(col_text):
                    data_cols.append([col_text, [c_ind_st, c_len], 'data'])
                if self.isDescription(col_text):
                    data_cols.append([col_text, [c_ind_st, c_len], 'description'])
            #print(data_cols)
            ind_description = 0
            if data_cols[0][2]=='description':
                ind_description = data_cols[0][1][0]
            ##
            data_titles = [None]*len(data_cols)
            for j in range(i-1, -1, -1):
                row = rows[j]
                cols = row.findAll('td')
                ## assert same number of columns
                n_col_j = 0
                for col in cols:
                    c_ind_st = n_col_j
                    c_len = 1
                    if col.has_attr('colspan'):
                        c_len = int(col['colspan'])
                    n_col_j += c_len
                if n_col_j!=n_col:
                    continue
                ##
                n_col_j = 0
                for col in cols:
                    c_ind_st = n_col_j
                    c_len = 1
                    if col.has_attr('colspan'):
                        c_len = int(col['colspan'])
                    n_col_j += c_len
                    col_text = col.text.strip()
                    if self.isColTitle(col_text) or (self.isDescription(col_text) and c_ind_st>ind_description):
                        #print(col_text)
                        for k in range(len(data_cols)):
                            if data_cols[k][2]=='data':
                                c_ind_st_k, c_len_k = data_cols[k][1]
                                if c_ind_st<=c_ind_st_k and c_ind_st+c_len>=c_ind_st_k+c_len_k:
                                    if data_titles[k]==None:
                                        data_titles[k] = col_text
                                    else:
                                        data_titles[k] = col_text + ' ' + data_titles[k]
            #print(data_titles)
            financialData.append([data_cols, data_titles])
        return financialData
        
        
    def removeParanthesis(self, title):
        if '(' in title and ')' in title:
            ind_0 = title.find('(')
            ind_1 = title.find(')')
            assert ind_0>=0 and ind_1>=0
            if ind_0<ind_1:
                title = title[:ind_0] + title[ind_1+1:]
        return title.strip()
        
    def replaceWord(self, title, word, word_replace):
        while word in title:
            ind_0 = title.find(word)
            ind_1 = ind_0 + len(word)
            title = title[:ind_0] + word_replace + title[ind_1:]
        return title
        
    def isTotalAssets(self, _title_):
        """ 0: No
            5: highest possibility
        """
        title = _title_.lower().strip()
        if title=='total assets':
            return 5
        elif 'total assets' in title:
            title = self.removeParanthesis(title)
            if title=='total assets':
                return 4
            else:
                self.logToFile('isTotalAssets: {} - not recognized.'.format(_title_))
        return 0        


    def isNetProfitOrLoss(self, _title_):
        """ 0: No
            5: highest possibility
        """
        title = _title_.lower().strip()
        title = self.replaceWord(title, '(loss)', 'loss')
        title = self.replaceWord(title, 'income', 'profit')
        title = title.replace(r' / ', ' ').replace(r' /', ' ').replace(r'/ ', ' ').replace(r'/', ' ')
        if title in ['net profit', 'net loss', 'net loss profit', 'net profit loss']:
            return 5
        elif 'net profit' in title:
            title = self.removeParanthesis(title)
            if title == 'net profit':
                return 4
            else:
                self.logToFile('isNetProfitOrLoss: {} - not recognized.'.format(_title_))
        elif 'net loss' in title:
            title = self.removeParanthesis(title)
            if title == 'net loss':
                return 4
            else:
                self.logToFile('isNetProfitOrLoss: {} - not recognized.'.format(_title_))
        
        return 0        
                        

    def isTotalLiabilities(self, _title_):
        """ 0: No
            5: highest possibility
        """
        title = _title_.lower().strip()
        if title=='total liabilities':
            return 5
        elif 'total liabilities' in title:
            title = self.removeParanthesis(title)
            if title=='total liabilities':
                return 4
            elif 'shareholder' in title or 'stockholder' in title:
                return 0
            else:
                self.logToFile('isTotalLiabilities: {} - not recognized.'.format(_title_))
        return 0        
        
    
    def isNetProfitPerShare(self, _title_):
        """ 0: No
            5: highest possibility
        """
        title = _title_.lower().strip()
        title = self.replaceWord(title, '(loss)', 'loss')
        title = self.replaceWord(title, 'income', 'profit')
        title = title.replace(r' / ', ' ').replace(r' /', ' ').replace(r'/ ', ' ').replace(r'/', ' ')
        title = self.removeParanthesis(title)
        
        if 'net loss profit per common share' in title or\
           'net profit per common share' in title or\
           'net loss per common share' in title or\
           'net profit loss per common share' in title:
            return 5
        else:
            self.logToFile('isNetProfitPerShare: {} - not recognized.'.format(_title_))
        return 0        


    def isNetProfitPerShareBasic(self, _title_):
        """ 0: No
            5: highest possibility
        """
        title = _title_.lower().strip()
        likelihood = self.isNetProfitPerShare(_title_)
        if likelihood>0:
            if 'dilut' not in title:
                return likelihood
        return 0        

    def isNetProfitPerShareDiluted(self, _title_):
        """ 0: No
            5: highest possibility
        """
        title = _title_.lower().strip()
        likelihood = self.isNetProfitPerShare(_title_)
        if likelihood>0:
            if 'dilut' in title:
                return likelihood
        return 0        



    def GetOrderAndYearlyForm(self, _title_):
        ##TODO: treat quarter labels as well
        order, yearlyform = None, None
        title = _title_.lower().strip()
        for year in self.yeartitles:
            if year in title:
                if 'three month' in title or '3 month' in title or 'quarter ended' in title:
                    yearlyform = self.ColLabels.M3Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
                elif 'nine month' in title or '9 month' in title:
                    yearlyform = self.ColLabels.M9Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
                elif 'twelve month' in title or '12 month' in title or 'year ended' in title:
                    yearlyform = self.ColLabels.M12Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
                elif 'six month' in title or '6 month' in title:
                    yearlyform = self.ColLabels.M6Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
                elif 'thirteen weeks' in title or '13 weeks' in title or\
                   'fourteen weeks' in title or '14 weeks' in title or\
                   'twelve weeks' in title or '12 weeks' in title:
                    yearlyform = self.ColLabels.W13Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
                elif 'thirty-nine weeks' in title or 'thirty nine weeks' in title or '39 weeks' in title or\
                   'thirty-eight weeks' in title or 'thirty eight weeks' in title or '38 weeks' in title or\
                   'fourty weeks' in title or '40 weeks' in title:
                    yearlyform = self.ColLabels.W39Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
                elif 'first quarter' in title or 'q1' in title or '1q' in title:
                    yearlyform = self.ColLabels.Q1Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        self.logToFile('GetOrderAndYearlyForm: data not recognized: \n{}'.format(_title_))
                    return [order, yearlyform]
                elif 'second quarter' in title or 'q2' in title or '2q' in title:
                    yearlyform = self.ColLabels.Q2Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        self.logToFile('GetOrderAndYearlyForm: data not recognized: \n{}'.format(_title_))
                    return [order, yearlyform]
                elif 'third quarter' in title or 'q3' in title or '3q' in title:
                    yearlyform = self.ColLabels.Q3Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        self.logToFile('GetOrderAndYearlyForm: data not recognized: \n{}'.format(_title_))
                    return [order, yearlyform]
                elif 'fourth quarter' in title or 'q4' in title or '4q' in title:
                    yearlyform = self.ColLabels.Q4Y
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        self.logToFile('GetOrderAndYearlyForm: data not recognized: \n{}'.format(_title_))
                    return [order, yearlyform]
                else:
                    yearlyform = self.ColLabels.YY
                    order_m = self.GetMonthlyOrder(title)
                    order = self.Y_orders[year]*12
                    if order_m!=None:
                        order += order_m
                    return [order, yearlyform]
        return None
                
    def GetMonthlyOrder(self, title):
        order_m = None
        if 'january' in title:
            order_m = 0
        elif 'february' in title:
            order_m = 1
        elif 'march' in title:
            order_m = 2
        elif 'april' in title:
            order_m = 3
        elif 'june' in title:
            order_m = 5
        elif 'july' in title:
            order_m = 6
        elif 'august' in title:
            order_m = 7
        elif 'september' in title:
            order_m = 8
        elif 'october' in title:
            order_m = 9
        elif 'november' in title:
            order_m = 10
        elif 'december' in title:
            order_m = 11
        elif 'may' in title:
            order_m = 4
        return order_m
        
        
    def ProcessFinancialDataStrings(self, financialData, lookFor):
        """ lookFor = self.descLabels.totalAssets etc.
        """
        self.descLabelsChecker = {self.descLabels.totalAssets:self.isTotalAssets, 
                                  self.descLabels.totalLiabilities:self.isTotalLiabilities,
                                  self.descLabels.netProfit:self.isNetProfitOrLoss,
                                  self.descLabels.netProfitPerShareBasic:self.isNetProfitPerShareBasic,
                                  self.descLabels.netProfitPerShareDiluted:self.isNetProfitPerShareDiluted}
                                  
        financialDataProcessed = []
        #print('financialData: ', financialData)
        if financialData==None:
            return financialDataProcessed
        for data_cols_titles in financialData:
            if data_cols_titles==None:
                continue
            data_cols, data_titles = data_cols_titles
            n_col = len(data_cols)
            assert len(data_titles)==n_col
            for i in range(n_col):
                if data_cols[i][2]=='description':
                    likelihood = self.descLabelsChecker[lookFor](data_cols[i][0])
                    if likelihood>0:
                        for j in range(n_col):
                            if data_cols[j][2]=='data':
                                if data_titles[j]==None:
                                    continue
                                order_form = self.GetOrderAndYearlyForm(data_titles[j])
                                if order_form!=None:
                                    order, yearlyform = order_form
                                    dataDic = {}
                                    dataDic['title'] = lookFor
                                    dataDic['likelihood'] = likelihood
                                    dataDic['data'] = self.numericStrToFloat(data_cols[j][0])
                                    dataDic['yearlyform'] = yearlyform
                                    dataDic['order'] = order
                                    financialDataProcessed.append(dataDic)
                
        return financialDataProcessed
        
    
    def ProcessFinancialDataChange(self, tables):
        total_change_list = []
        for i in range(len(tables)):
            table = tables[i]

            asset_data = self.HtmlTable_FindFinancialTermData(table, 'total assets')
            liability_data = self.HtmlTable_FindFinancialTermData(table, 'total liabilities')
            profit_data = self.HtmlTable_FindFinancialTermData(table, 'net profit')
            income_data = self.HtmlTable_FindFinancialTermData(table, 'net income')
            loss0_data = self.HtmlTable_FindFinancialTermData(table, 'net loss')
            loss1_data = self.HtmlTable_FindFinancialTermData(table, 'net (loss)')

            asset_data_processed = self.ProcessFinancialDataStrings(asset_data, lookFor=self.descLabels.totalAssets)
            liability_data_processed = self.ProcessFinancialDataStrings(liability_data, lookFor=self.descLabels.totalLiabilities)
            profit_data_processed = self.ProcessFinancialDataStrings(profit_data, lookFor=self.descLabels.netProfit)
            profit_data_processed.extend(self.ProcessFinancialDataStrings(income_data, lookFor=self.descLabels.netProfit))
            profit_data_processed.extend(self.ProcessFinancialDataStrings(loss0_data, lookFor=self.descLabels.netProfit))
            profit_data_processed.extend(self.ProcessFinancialDataStrings(loss1_data, lookFor=self.descLabels.netProfit))
            
            
            asset_change = self.CalculatePercentChange(asset_data_processed, 
                w_pp=self.weightingFactors_pp['asset'], w_pm=self.weightingFactors_pm['asset'], 
                w_mp=self.weightingFactors_mp['asset'], w_mm=self.weightingFactors_mm['asset'])
            liabilitiy_change = self.CalculatePercentChange(liability_data_processed, 
                w_pp=self.weightingFactors_pp['liability'], w_pm=self.weightingFactors_pm['liability'], 
                w_mp=self.weightingFactors_mp['liability'], w_mm=self.weightingFactors_mm['liability'])
            profit_change = self.CalculatePercentChange(profit_data_processed,
                w_pp=self.weightingFactors_pp['profit'], w_pm=self.weightingFactors_pm['profit'], 
                w_mp=self.weightingFactors_mp['profit'], w_mm=self.weightingFactors_mm['profit'])
            
            if self.vbose:
                print(asset_data, '\n', '-'*30)
                print(liability_data, '\n', '-'*30)
                print(profit_data, '\n', '-'*30)
                print(income_data, '\n', '-'*30)
                print(loss0_data, '\n', '-'*30)
                print(loss1_data, '\n', '-'*30)
                print('='*80)
                print(asset_data_processed, '\n', 'asset_change: ', asset_change, '\n', '-'*30)
                print(liability_data_processed, '\n', 'liabilitiy_change: ', liabilitiy_change, '\n', '-'*30)
                print(profit_data_processed, '\n', 'profit_change: ', profit_change, '\n', '-'*30)
                print('=<>='*30)
        
            n_elem = 0
            total_change = None
            if asset_change!=None:
                if total_change==None:
                    total_change = 0.0
                total_change += asset_change
                n_elem += 1
            if liabilitiy_change!=None:
                if total_change==None:
                    total_change = 0.0
                total_change += liabilitiy_change
                n_elem += 1
            if profit_change!=None:
                if total_change==None:
                    total_change = 0.0
                total_change += profit_change
                n_elem += 1
            
            if total_change!=None and n_elem>0:
                total_change /= n_elem
            total_change_list.append(total_change)
            
        return np.array(total_change_list, dtype=np.float)
        
            

    def CalculatePercentChange(self, data_processed, w_pp=None, w_pm=None, w_mp=None, w_mm=None):
        if data_processed==None or len(data_processed)==0:
            return None
            
        yearlyforms = {}
        
        for i in range(len(data_processed)):
            yf_i = data_processed[i]['yearlyform']
            if yf_i not in yearlyforms:
                yearlyforms[yf_i] = {data_processed[i]['order']:[(data_processed[i]['likelihood'], data_processed[i]['data'])]} 
            else:
                if data_processed[i]['order'] in yearlyforms[yf_i]:
                    yearlyforms[yf_i][data_processed[i]['order']].append((data_processed[i]['likelihood'], data_processed[i]['data']))
                else:
                    yearlyforms[yf_i][data_processed[i]['order']] = [(data_processed[i]['likelihood'], data_processed[i]['data'])]
                    
        ##keep maximum likelihood and drop redundencies
        for yf in yearlyforms:
            for order in yearlyforms[yf]:
                yearlyforms[yf][order] = list(set(yearlyforms[yf][order]))
                if len(yearlyforms[yf][order])>1:
                    lk_max = 0
                    likelihood_data_list = yearlyforms[yf][order]
                    for lk, d in likelihood_data_list:
                        if lk_max<lk:
                            lk_max = lk
                    for j in range(len(likelihood_data_list)-1, -1, -1):
                        if likelihood_data_list[j][0]<lk_max:
                            del likelihood_data_list[j]
                if len(yearlyforms[yf][order])>1:
                    self.logToFile('CalculatePercentChange: identical orders:\n {}.'.format(data_processed))
                    return None
                
        ##calculate percent change
        p_ch = 0
        for yf in yearlyforms:
            yf_unordered = [(order, yearlyforms[yf][order]) for order in yearlyforms[yf]]
            yf_ordered = sorted(yf_unordered, key=lambda x: x[0])
            #print('yf_ordered: ', yf_ordered)
            yf_ordered = [yf_ordered[i][1][0][1] for i in range(len(yf_ordered))]
            for i in range(len(yf_ordered)-1):
                p_ch_i = 0.0
                if abs(yf_ordered[i])>=1.0:
                    p_ch_i = (yf_ordered[i+1]-yf_ordered[i])/abs(yf_ordered[i])
                else:
                    p_ch_i = yf_ordered[i+1]-yf_ordered[i]

                if yf_ordered[i]>=0.0 and p_ch_i>=0.0 and w_pp!=None:
                    p_ch_i *= w_pp
                elif yf_ordered[i]>=0.0 and p_ch_i<0.0 and w_pm!=None:
                    p_ch_i *= w_pm
                elif yf_ordered[i]<0.0 and p_ch_i>=0.0 and w_mp!=None:
                    p_ch_i *= w_mp
                elif yf_ordered[i]<0.0 and p_ch_i<0.0 and w_mm!=None:
                    p_ch_i *= w_mm

                p_ch += p_ch_i
        
        return p_ch


    def SetLastDateChecked(self, last_checked):
        self.datetime_last_checked = last_checked

    def CheckSECReportsUntilLastDateChecked(self, printall=True, printpositive=True, alarmOnPositive=True):
        page_no=0
        stop = False
        last_checked = None
        while not stop:
            data_list = self.GetRecentFilings_DirectWebScrap(page_no=page_no, displayTable=False)

            ind_hasItem = [i for i in range(len(data_list)) if data_list[i]['hasItems']==True] 
            #print('Indices with requested item :', ind_hasItem)
            if page_no==0:                
                last_checked = self.GetAcceptedDatetime(data_list[0]['accepted'])

            for item_no in range(len(ind_hasItem)):
                if self.GetAcceptedDatetime(data_list[ind_hasItem[item_no]]['accepted'])<=self.datetime_last_checked:
                    stop = True
                    break
                    
                self.logToFile("{}\n{}\n{}".format('='*80, data_list[ind_hasItem[item_no]]['name'], data_list[ind_hasItem[item_no]]['link']))

                txtFile = self.GetItemTextFile_WebScrap(data_list, ind_hasItem[item_no])

                tables = self.TxtFile_ExtractTables(txtFile, displayTables=False)

                total_changes = self.ProcessFinancialDataChange(tables)
                total_changes = np.ma.masked_invalid(total_changes)

                if printall:
                    print(total_changes)
                    print(data_list[ind_hasItem[item_no]]['name'])
                    print(data_list[ind_hasItem[item_no]]['accepted'])
                    print(data_list[ind_hasItem[item_no]]['link'])
                    print('Mean: ', total_changes.mean())
                    print('-'*50)
                if printpositive:
                    if isinstance(total_changes.mean(), float) and total_changes.mean()>0.0:
                        print(total_changes)
                        print(data_list[ind_hasItem[item_no]]['name'])
                        print(data_list[ind_hasItem[item_no]]['accepted'])
                        print(data_list[ind_hasItem[item_no]]['link'])
                        print('Mean: ', total_changes.mean())
                        print('='*50)
                        
                        if alarmOnPositive:
                            self.SoundAlarm()
                        
            page_no += 1
        return last_checked


    def CheckSECEvery(self, minute=2, second=0, printall=True, printpositive=True, stopDatetime=None):
        dt = int(60*minute + second)
            
        while True:
            last_checked = self.CheckSECReportsUntilLastDateChecked(printall=printall, printpositive=printpositive)
            
            if last_checked!=None:
                self.SetLastDateChecked(last_checked)
                
            time.sleep(dt)
            
            if stopDatetime!=None:
                if datetime.datetime.now()>stopDatetime:
                    break
            

    def SoundAlarm(self, n=1):
        ##launches a wndow
        #os.system('ffplay {}'.format(self.mp3alarm_path))
        ##pure terminal
        for i in range(n):
            os.system('mpg123 {}'.format(self.mp3alarm_path))
            time.sleep(1)
    
        
    def FeedItem_GetUrl(self, item_no):
        url_dataPage = self.feed[ "items" ][item_no]['link']
        if self.vbose:
            print('link to data page: \n', url_dataPage)
        return url_dataPage
            
    def FeedItem_GetUrlHTMLPage(self, item_no):
        url_dataPage = self.FeedItem_GetUrl(item_no)
        html = urllib.request.urlopen(url_dataPage).read().decode('utf-8')
        return html
        
            
    def HTMLPage_GetExcelLink(self, html):
        if html.find('Interactive Data')>=0:
            #print('Excel file available.')
            bs = BeautifulSoup(html, 'lxml')
            #link_interactiveData = bs.find(lambda tag: tag.name=='a' and tag.has_attr('href') and tag.has_attr('id') and tag['id']=="interactiveDataBtn") 
            link_interactiveData = bs.find('a', href=True, text="\u00a0Interactive Data") 
            #print(link_interactiveData)
            link_interactiveData = self.SECLink + link_interactiveData['href']
            
            if self.vbose:
                print('link to interactive data page: \n', link_interactiveData)

            html_excelPage = urllib.request.urlopen(link_interactiveData).read().decode('utf-8')
            bs = BeautifulSoup(html_excelPage, 'lxml')
            link_excel = bs.find('a', href=True, text="View Excel Document") 
            #print(link_excel)
            link_excel = self.SECLink + link_excel['href']
            
            if self.vbose:
                print('link to excel file: \n', link_excel)
            
            return link_excel
        else:
            return None
            

    def HtmlPage_extractTxtFile(self, html, save_file=False, filepath='other/file.html'):    

        bs = BeautifulSoup(html, 'lxml')
        table = bs.find(lambda tag: tag.name=='table' and tag.has_attr('summary') and tag['summary']=="Document Format Files") 
        rows = table.findAll(lambda tag: tag.name=='tr')
        rows = str(rows)
        link_txtData_start = rows.find('Complete submission text file')
        link_txtData_end = rows.find('.txt</a></td>', link_txtData_start)
        link_txtData = rows[link_txtData_start:link_txtData_end]

        link_txtData_start = link_txtData.find('a href="') + 8
        link_txtData_end = link_txtData.find('.txt"', link_txtData_start) + 4
        link_txtData = 'https://www.sec.gov'+link_txtData[link_txtData_start:link_txtData_end]
        if self.vbose:
            print('link to complete txt submission: \n', link_txtData)

        file_txt = urllib.request.urlopen(link_txtData).read().decode('utf-8')

        if save_file:
            with open(filepath, 'w') as file_:
                file_.write(file_txt)
        return file_txt


    def TxtFile_ExtractTables(self, file_txt, filterData=True, displayTables=True):
        bs = BeautifulSoup(file_txt, 'lxml')
        tables = bs.findAll(lambda tag: tag.name=='table') 

        if filterData:
            look_for = ['current assets', 'net loss', 'net income', 'total assets', 'total liabilities', 'net profit', \
                        'total debt', 'net (loss)']
            tables_filtered = []
            for table in tables:
                for financial_term in look_for:
                    table_lower = str(table).lower()
                    if table_lower.find(financial_term)>=0:
                        tables_filtered.append(table)
                        if displayTables:
                            display(HTML(str(table)))
            return tables_filtered

        else:
            return tables
        
        


    def HtmlTable_FindFinancialTermData(self, table, financial_term='total assets'):
        return self.HtmlTable_GetTermAndPeriod(table, financial_term.lower())
                
            
    
    def GetCompanyReportsFeed(self, CIK, page_no=0):
        ## CIK = '0001338042'

        count = self.count
        start = page_no*count

        edgar_cik_rss_url = 'https://www.sec.gov/cgi-bin/browse-edgar?action=getcompany&CIK={}&CIK={}&type=&dateb=&owner=include&start={}&count={}&output=atom'.format(CIK, CIK, start, count)
        feed = feedparser.parse( edgar_cik_rss_url )
        
        return feed
    
    ##direct web scraping (updated realtime as opposed to rss feeds updated each 10 minutes)
    def GetTableHeaders(self, table):
        heads = table.findAll(lambda tag: tag.name=='th')
        headers =[]
        for h in heads:
            headers.append(h.text.strip().lower())
        return headers
            
    def IdentifyRecentReportedTable(self, tables):
        table_good = None
        for table in tables:
            if table==None:
                continue
            headers = self.GetTableHeaders(table)
            if headers==None:
                continue
            headers_good = ['file', 'formats', 'description', 'accepted', 'filing date', r'file/film no'] 
            n_good =  0
            for h in headers_good:
                if h in headers:
                    n_good += 1
            #print(n_good)

            if n_good>3:
                table_good = table
                break
        return table_good

    def GetRecentFilings_DirectWebScrap(self, page_no=0, displayTable=False):
        count = self.count
        form_type = self.form_type
        form_items = self.form_items

        start = page_no*count

        url_recent = 'https://www.sec.gov/cgi-bin/browse-edgar?company=&CIK=&type={}&owner=include&count={}&action=getcurrent'.format(form_type, count)
        if page_no>0:
            url_recent = 'https://www.sec.gov/cgi-bin/browse-edgar?action=getcurrent&datea=&dateb=&company=&type={}&SIC=&State=&Country=&CIK=&owner=include&accno=&start={}&count={}'.format(form_type, start, count)

        html = urllib.request.urlopen(url_recent).read().decode('utf-8')

        bs = BeautifulSoup(html, 'lxml')
        tables = bs.findAll(lambda tag: tag.name=='table')# and tag.has_attr('summary') and tag['summary']=="Document Format Files") 
    

        table = self.IdentifyRecentReportedTable(tables)

        headers = self.GetTableHeaders(table)
        headers_ind = dict(zip(headers, range(len(headers))))
        
        if self.vbose:
            print('headers: ', headers)
            print('headers_ind: ', headers_ind)
        
        if displayTable:
            display(HTML(str(table)))
        rows = table.findAll(lambda tag: tag.name=='tr')

        #print('No rows: ', len(rows))
        assert len(rows)==2*count+1

        data_list = []
        for i in range(count):
            row = rows[2*i+1]
            cols = row.findAll('td')
            name = cols[headers_ind['description']].text
            
            link = 'https://www.sec.gov'+cols[headers_ind['description']].find('a', href=True)['href']

            row = rows[2*i+2]
            cols = row.findAll('td')
            description = cols[headers_ind['description']].text
            hasItems = False
            if form_items==None or len(form_items)==0:
                hasItems = True
            else:
                for f_item in form_items:
                    if f_item in description:
                        hasItems = True
            
            txt_url = 'https://www.sec.gov'+cols[headers_ind['formats']].find('a', href=True, text="[text]")['href']

            accepted = cols[headers_ind['accepted']].text
            
            data_list.append({'name':name, 'link':link, 'description':description,\
                              'accepted':accepted, 'txt_url':txt_url, 'hasItems':hasItems})
            
            if self.vbose:
                print('{}: {}\n{}\n{}{}\n{}\n{}\n'.format(i, name, link, description, accepted, txt_url, hasItems))
        
        self.WebScrap_RecentForms_DataList = data_list
        return data_list

    
    def GetAcceptedDatetime(self, accepted_datetime_str):
        return datetime.datetime.strptime(accepted_datetime_str, '%Y-%m-%d%H:%M:%S')
    
            
    def GetItemTextFile_WebScrap(self, data_list, item_no, save_file=False, filepath='other/file.html'):
        txt_url = data_list[item_no]['txt_url']
        file_txt = urllib.request.urlopen(txt_url).read().decode('utf-8')

        if save_file:
            with open(filepath, 'w') as file_:
                file_.write(file_txt)
        return file_txt
            
            
            
            
            
                    
