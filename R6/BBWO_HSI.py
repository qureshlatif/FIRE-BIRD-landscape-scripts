import scipy
from scipy import spatial
from scipy import stats
from numpy import *
from scipy import special

class Mahalonobis_HSI():
    def mahalanobis_sqdist(self, d, mean, Sigma):
        Sigma_inv = linalg.cholesky(Sigma)
        xdiff = d - mean
        d2 = sum(linalg.solve(Sigma_inv, xdiff.T)**2, axis = 0)
        return d2

    def Mahal(self, Cal,Lnd):
        mns = mean(Cal,0)
        sds = std(Cal,0,ddof = 1)
        cal = (Cal - mns)/sds
        lnd = (Lnd - mns)/sds
        cnt = mean(cal, 0)
        cv = cov(cal.T)
        dgf = size(cal,1)
        D_lnd = self.mahalanobis_sqdist(lnd, cnt, cv)
        p_lnd =  1 - special.chdtr(dgf, D_lnd)
        return p_lnd

    def Mxnt3v_scr(self, cosasp, dr,loccc40): 
      #Define features and parameters
      l_cosasp =array([0.9381002574647688, -1, 1])
      
      l_dnbr = array([10.786746325527512, -590.47619629, 1024.12268066])
      dnbr = copy(dr)
      dnbr[dnbr > l_dnbr[2]] = l_dnbr[2]
      dnbr[dnbr < l_dnbr[1]] = l_dnbr[1]
      
      l_loccc40 = array([1.1312580836869854, 0.0, 1.0])
      
      l_cosasp2 = array([-0.3980392534351936, 0.0, 1.0])
      cosasp2 = copy(cosasp)**2
      
      l_dnbr2 = array([-6.415155660077871, 0.0, 1048827.2650422244])
      dnbr2 = copy(dnbr)**2
      dnbr2[dnbr2 > l_dnbr2[2]] = l_dnbr2[2]
      dnbr2[dnbr2 < l_dnbr2[1]] = l_dnbr2[1]
      
      l_loccc402 = array([-1.9494267010037454, 0.0, 1.0])
      loccc402 = copy(loccc40)**2
      
      l_ca_x_lcc40 = array([-0.08235355467277128, -1.0, 1.0])
      ca_x_lcc40 = copy(cosasp) * copy(loccc40)
      ca_x_lcc40[ca_x_lcc40 > l_ca_x_lcc40[2]] = l_ca_x_lcc40[2]
      ca_x_lcc40[ca_x_lcc40 < l_ca_x_lcc40[1]] = l_ca_x_lcc40[1]
      
      l_dnbr_x_lcc40 = array([3.6550086289032, -213.7978722982045, 1024.0771640968703])
      dnbr_x_lcc40 = copy(dnbr) * copy(loccc40)
      dnbr_x_lcc40[dnbr_x_lcc40 > l_dnbr_x_lcc40[2]] = l_dnbr_x_lcc40[2]
      dnbr_x_lcc40[dnbr_x_lcc40 < l_dnbr_x_lcc40[1]] = l_dnbr_x_lcc40[1]
      
      linPN = 8.12466828053413
      densNorm = 2846.1535672273485
      entropy = 8.807545807433492
      
      exponent = (l_cosasp[0] * ((cosasp - l_cosasp[1]) /( l_cosasp[2] - l_cosasp[1])) +
                     l_dnbr[0] * ((dnbr - l_dnbr[1]) / (l_dnbr[2] - l_dnbr[1])) +
                     l_loccc40[0] * ((loccc40 - l_loccc40[1]) / (l_loccc40[2] - l_loccc40[1])) +
                     l_cosasp2[0] * ((cosasp2 - l_cosasp2[1]) / (l_cosasp2[2] - l_cosasp2[1])) +
                     l_dnbr2[0] * ((dnbr2 - l_dnbr2[1]) / (l_dnbr2[2] - l_dnbr2[1])) +
                     l_loccc402[0] * ((loccc402 - l_loccc402[1]) / (l_loccc402[2] - l_loccc402[1])) +
                     l_ca_x_lcc40[0] * ((ca_x_lcc40 - l_ca_x_lcc40[1]) / (l_ca_x_lcc40[2] - l_ca_x_lcc40[1])) +
                     l_dnbr_x_lcc40[0] * ((dnbr_x_lcc40 - l_dnbr_x_lcc40[1]) / (l_dnbr_x_lcc40[2] - l_dnbr_x_lcc40[1]))
                     ) - linPN
      
      mx_raw  = exp(exponent)/densNorm
      HSI = (mx_raw*exp(entropy))/(1+mx_raw*exp(entropy))
      return HSI
    
    def Mxntbrn_scr(self, dr):
        l_dnbr = array([5.393373689029553, -590.47619629, 1024.12268066])
        dnbr = copy(dr)
        dnbr[dnbr>l_dnbr[2]] = l_dnbr[2]
        dnbr[dnbr<l_dnbr[1]] = l_dnbr[1]
        linPN = 5.393373689029553
        densNorm = 1369.5883369712144
        entropy = 8.899831912829168
        exponent = l_dnbr[0]*((dnbr-l_dnbr[1])/(l_dnbr[2]-l_dnbr[1])) - linPN
        mx_raw  = exp(exponent)/densNorm
        HSI = (mx_raw*exp(entropy))/(1+mx_raw*exp(entropy))
        return HSI


    def expit(self, x):
        return exp(x)/(1+exp(x))

    def HSI_calc(self, cosasp,dnbr,loccc40,lndcc40):
      SG_wlr = self.expit(-3.130234119 + 0.003704316*dnbr + 3.154616050*lndcc40)
      TP_wlr = self.expit(-2.485159177 + 0.006151492*dnbr)
      TB_wlr = self.expit(-0.743163 + 0.001546*dnbr)
      Maxent_3v = self.Mxnt3v_scr(cosasp,dnbr,loccc40)
      Maxent_brn = self.Mxntbrn_scr(dnbr)
      HSIs = array([SG_wlr,TP_wlr,TB_wlr,Maxent_3v,Maxent_brn])
      return HSIs 
  
