"""
This script retrieves the stage out parameters for your crab.cfg using the Phedex LFN2PFN API, if you want to stage out into group space /store/group.

You need to be in a CMSSW release environment (SL5 releases only!) before using this script.

For example:
    cmsrel CMSSW_x_x_x
    cd CMSSW_x_x_x
    cmsenv
    

Author: Manuel Giffels <giffels@physik.rwth-aachen.de> (2010)
"""

from xml.dom.minidom import parse

import urllib

try:
    import json
except:
    import simplejson as json
    
    
class PhEDExDatasvcInfo:

    def __init__(self,datasvc_url="http://cmsweb.cern.ch/phedex/datasvc",format="json",instance="prod"):
        self.format=format
        self.url = datasvc_url+'/'+self.format+'/'+instance

    def __GetDomFromPhedex(self,parameters,api):  

        SupportedAPIs = ['lfn2pfn']

        if api not in SupportedAPIs:
            msg="API %s is currently not supported." % (api)
            raise NotImplementedError,msg
    
        parameters = urllib.urlencode(parameters)

        try:
            urlresults = urllib.urlopen(self.url+'/'+api+'?'+parameters)
            urlresults = parse(urlresults)
        except IOError:
            msg="Unable to access PhEDEx Data Service at %s with given API %s" % (self.url,api)
            raise IOError,msg

        return urlresults

    def __GetJSONFromPhedex(self,parameters,api):

        SupportedAPIs = ['lfn2pfn']

        if api not in SupportedAPIs:
            msg="API %s is currently not supported." % (api)
            raise NotImplementedError,msg
    
        parameters = urllib.urlencode(parameters)

        try:
            urlresults = urllib.urlopen(self.url+'/'+api+'?'+parameters)
            urlresults = json.load(urlresults)
        except IOError:
            msg="Unable to access PhEDEx Data Service at %s with given API %s" % (self.url,api)
            raise IOError,msg

        return urlresults

    def GetPFNFromLFN(self,node,lfn,protocol='srmv2'):
        parameters = {'node' : node , 'lfn': lfn , 'protocol': protocol}

        if self.format == "xml":

            domresults = self.__GetDomFromPhedex(parameters,'lfn2pfn')
            
            result = domresults.getElementsByTagName('phedex')

            if not result:
                return []

            result = result[0]
            pfn = None
            mapping = result.getElementsByTagName('mapping')

            for i in mapping:
                pfn=i.getAttribute("pfn")
                if pfn:
                    return pfn

        elif self.format == "json":
            jsondict = self.__GetJSONFromPhedex(parameters,'lfn2pfn')

            if not jsondict.has_key('phedex') or not jsondict['phedex'].has_key('mapping'):
                print jsondict
                return None
            
            for i in jsondict['phedex']['mapping']:
                try:
                    return i['pfn']
                except:
                    return None
        else:
            return None
        
def main():
    import getopt,sys
    from urlparse import urlparse

    usage = "python GetStorageCfgFromTFC.py --site=T2_XY_ABCDE --lfn=/store/user/AliBaba/..."

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help", "site=","lfn="])
    except getopt.GetoptError, err:
        print str(err) 
        print usage
        sys.exit(-1)

    cmsname = None
    lfn=None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print usage
            sys.exit(-1)
        elif opt in ("--site=="):
            cmsname = arg.strip()
        elif opt in ("--lfn=="):
            lfn = arg.strip()
        else:
            assert False, "unhandled option"

    if cmsname==None or lfn==None:
        print usage
        sys.exit(-1)

    phedex = PhEDExDatasvcInfo()
    
    srmurl = phedex.GetPFNFromLFN(node=cmsname,lfn=lfn)

    if srmurl:
        result = urlparse("http://"+srmurl.split('://')[1])#srm is not a valid schema, to be fixed
        print "storage_element="+result.hostname
        print "storage_port="+str(result.port)
        print "storage_path="+result.path+'?'+result.query.replace(lfn,"")
        print "user_remote_dir="+lfn
    else:
        print "Association not possible. Please, check input parameters."

if __name__ == '__main__':
    main()




