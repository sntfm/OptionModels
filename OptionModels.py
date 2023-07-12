import numpy as np
from scipy.stats import norm

class GeneralizedBlackScholes:
    pass

class Black76:
    pass

class BlackScholes:
    class OptionBS:
        def __init__(self, S= None, K= None, V= None, r= None, T= None, t=0, opt= 'call', value= None):
            # market input
            self.S, self.K, self.V, self.r, self.type = S, K, V, r, opt
            self.T, self.t = T/365, t/365
            self.tau = self.T - self.t

            # d1d2
            self.d1= (np.log(self.S/self.K) + (self.r + (self.V**2)/2)*self.tau)\
                        /(self.V*np.sqrt(self.tau))
            self.d2= self.d1 - self.V*np.sqrt(self.tau)

            # value and greeks
            self.value= value
            self.delta= None
            self.gamma= None
            self.vega= None
            self.theta= None
            self.rho= None

    def __init__(self, S= None, K= None, V= None, r= None, T= None, t= 0, opt= 'call', value= None):
            self.option= self.OptionBS(S= S, K= K, V= V, r= r, T= T, t= t,  opt= opt, value= value)
            self.value(self.option)
            self.delta(self.option)
            self.gamma(self.option)
            self.vega(self.option)
            self.theta(self.option)
            self.rho(self.option)

    def value(self, option):
        if option.type == 'call':    
            valueCall = option.S * norm.cdf(option.d1, 0, 1) - \
                option.K * np.exp(-option.r * option.tau) * norm.cdf(option.d2, 0, 1)
            option.value= valueCall
            return valueCall
        elif option.type == 'put':
            valuePut = option.K * np.exp(-option.r * option.tau) * norm.cdf(-option.d2, 0, 1) - \
                option.S * norm.cdf(-option.d1, 0, 1)
            option.value= valuePut
            return valuePut

    def delta(self, option):
        if option.type == 'call':
            deltaCall= norm.cdf(option.d1, 0, 1)
            option.delta= deltaCall
            return deltaCall
        elif option.type == 'put':
            deltaPut= norm.cdf(option.d1, 0, 1) - 1
            option.delta= deltaPut
            return deltaPut

    def gamma(self, option):
        gamma= norm.pdf(option.d1, 0, 1)/(option.S*option.V*np.sqrt(option.tau))
        option.gamma= gamma
        return gamma
    
    def vega(self, option):
        vega= option.S * norm.pdf(option.d1, 0, 1) * np.sqrt(option.tau)
        option.vega= vega/100
        return vega/100
    
    def theta(self, option):
        if option.type == 'call':
            thetaCall= -(option.S * norm.pdf(option.d1, 0, 1) * option.V)/(2*np.sqrt(option.tau)) \
                  - option.r * option.K * np.exp(-option.r * option.tau) * norm.cdf(option.d2, 0, 1)
            option.theta= thetaCall/365
            return thetaCall/365
        elif option.type == 'put':
            thetaPut= -(option.S * norm.pdf(option.d1, 0, 1) * option.V/(2*np.sqrt(option.tau))) \
                  + option.r * option.K * np.exp(-option.r * option.tau) * norm.cdf(-option.d2, 0, 1)
            option.theta= thetaPut/365
            return thetaPut/365

    def rho(self, option):
        if option.type == 'call':
            rhoCall= option.K * option.tau * np.exp(-option.r * option.tau) * norm.cdf(option.d2, 0, 1)
            option.rho= rhoCall/100
            return rhoCall/100
        elif option.type == 'put':
            rhoPut= -option.K * option.tau * np.exp(-option.r * option.tau) * norm.cdf(-option.d2, 0, 1)
            option.rho= rhoPut/100
            return rhoPut/100

class GarmanKohlhagen:
    class OptionGK:
        def __init__(self, S= None, K= None, V= None, rd= None, rf= None, T= None, t=0, opt= 'call', value= None):
            # market input
            self.S, self.K, self.V, self.rd, self.rf, self.type = S, K, V, rd, rf, opt
            self.T, self.t = T/365, t/365
            self.tau = self.T - self.t

            # d1d2
            self.d1 = (np.log(self.S/self.K) + (self.rd - self.rf + (self.V**2)/2)*self.tau)\
                        /(self.V*np.sqrt(self.tau))
            self.d2 = self.d1 - self.V*np.sqrt(self.tau)

            # value and greeks
            self.value= value
            self.delta= None
            self.gamma= None
            self.vega= None
            self.theta= None
            self.rho_d= None
            self.rho_f= None

    def __init__(self, S= None, K= None, V= None, rd= None, rf= None, T= None, t= 0, opt= 'call', value= None):
        self.option= self.OptionGK(S= S, K= K, V= V, rd= rd, rf= rf, T= T, t= t,  opt= opt, value= value)
        self.value(self.option)
        self.delta(self.option)
        self.gamma(self.option)
        self.vega(self.option)
        self.theta(self.option)
        self.rho_d(self.option)
        self.rho_f(self.option)

    def value(self, option):
        if option.type == 'call':    
            valueCall = option.S * np.exp(-option.rf*option.tau) * norm.cdf(option.d1, 0, 1)\
                - option.K * np.exp(-option.rd*option.tau) * norm.cdf(option.d2, 0, 1)
            option.value= valueCall
            return valueCall
        elif option.type == 'put':
            valuePut = option.K * np.exp(-option.rd*option.tau) * norm.cdf(-option.d2, 0, 1)\
                - option.S * np.exp(-option.rf*option.tau) * norm.cdf(-option.d1, 0, 1)
            option.value= valuePut
        return valuePut

    def delta(self, option):
        if option.type == 'call':
            deltaCall= np.exp(-option.rf*option.tau) * norm.cdf(option.d1, 0, 1)
            option.delta= deltaCall
            return deltaCall
        elif option.type == 'put':
            deltaPut= np.exp(-option.rf*option.tau) * norm.cdf(option.d1, 0, 1) - 1
            option.delta= deltaPut
            return deltaPut

    def gamma(self, option):
        gamma= (np.exp(-option.rf*option.tau) * norm.pdf(option.d1, 0, 1))\
            /(option.S*option.V*np.sqrt(option.tau))
        option.gamma= gamma
        return gamma   

    def vega(self, option):
        vega= np.exp(-option.rf*option.tau) * option.S * norm.pdf(option.d1, 0, 1) * np.sqrt(option.tau)
        option.vega= vega/100
        return vega/100

    def theta(self, option):
        if option.type == 'call':
            thetaCall= option.S * option.rf * np.exp(-option.rf*option.tau) * norm.cdf(option.d1, 0, 1)\
                - option.K * option.rd * np.exp(-option.rd*option.tau) * norm.cdf(option.d2, 0, 1)\
                - (option.S * option.V * np.exp(-option.rf*option.tau) * norm.pdf(option.d1, 0, 1))/(2*np.sqrt(option.tau))
            option.theta= thetaCall/365
            return thetaCall/365
        elif option.type == 'put':
            thetaPut= - option.S * option.rf * np.exp(-option.rf*option.tau) * norm.cdf(-option.d1, 0, 1)\
                + option.K * option.rd * np.exp(-option.rd*option.tau) * norm.cdf(-option.d2, 0, 1)\
                - (option.S * option.V * np.exp(-option.rf*option.tau) * norm.pdf(option.d1, 0, 1))/(2*np.sqrt(option.tau))
            option.theta= thetaPut/365
            return thetaPut/365

    def rho_d(self, option):
        if option.type == 'call':
            rhoCall= option.K * option.tau * np.exp(-option.rd * option.tau) * norm.cdf(option.d2, 0, 1)
            option.rho_d= rhoCall/100
            return rhoCall/100
        elif option.type == 'put':
            rhoPut= - option.K * option.tau * np.exp(-option.rd * option.tau) * norm.cdf(-option.d2, 0, 1)
            option.rho_d= rhoPut/100
            return rhoPut/100

    def rho_f(self, option):
        if option.type == 'call':
            rhoCall= - option.S * option.tau * np.exp(-option.rf * option.tau) * norm.cdf(option.d1, 0, 1)
            option.rho_f= rhoCall/100
            return rhoCall/100
        elif option.type == 'put':
            rhoPut= option.S * option.tau * np.exp(-option.rf * option.tau) * norm.cdf(-option.d1, 0, 1)
            option.rho_f= rhoPut/100
            return rhoPut/100 

class Heston:
    pass

class CoxRossRubenstein:
    pass

class ImpliedVolatility:
    pass
    # def impvol(self, S, K, r, T, opt, value, t=0, maxiter= 50, tol= 0.01, initvol= 0.5, alpha= 0.01):
    #     'Newton Raphson fast approximation'
    #     T= T/365
    #     t= t/365
    #     impvol= 0.19
    #     i= 0

    # while i < maxiter:
    #     option= self.Option(S= S, K= K, r= r, T= T, V= impvol, opt= opt)
    #     self.value(option)
    #     print(option.value)
    #     diff= -value + self.value(option)
    #     vega= self.vega(option)
    #     print(i, diff, vega, impvol, self.value(option))

    #     if abs(diff) <= tol:
    #         break
    #     print(diff/(100*vega))
    #     impvol += alpha*diff/(vega)
    #     i += 1
        
    # return impvol


    # def find_iv_newton(c_p, S, K, r, t, market_price):
    #     _sigma = 0.5
    #     for i in range(MAX_TRY):
            # _bs_price = bs_price(c_p,S, K, r, t, sigma=_sigma)
            # diff = market_price - _bs_price
        #     # vega = S*N_prime(d1)*sqrt(t)
        #     if abs(diff) < ONE_CENT:
        #         return _sigma
        #     _sigma += diff/vega
        # return _sigma

# BS smoke test
if __name__ == '__main__' and False:
    spot=90.5
    strike= 92.0
    vol= 0.19
    r= 0.053
    T= 85
    opt= 'call'

    bs= BlackScholes(S= spot, K= strike, V= vol, r= r, T= T, opt= opt)
    print(f'''
    Reference:
    Spot:   {bs.option.S}    Strike: {bs.option.K}
    Vol:    {bs.option.V}    Tau:{round(bs.option.tau,4)}
    R:      {bs.option.r}         

    <Option>  

    Type:  {bs.option.type.upper()}
    Value: {round(bs.option.value, 3)}

    <Greeks>

    Delta:      {round(bs.option.delta,3)}
    Gamma:      {round(bs.option.gamma,3)}
    Vega:       {round(bs.option.vega,3)}
    Theta:      {round(bs.option.theta,3)}
    Rho:        {round(bs.option.rho,3)}
    ''')

# GK smoke test
if __name__ == '__main__' and False:
    spot= 90.5
    strike= 92.0
    vol= 0.19
    rd= 0.053
    rf= 0.025
    T= 85
    opt= 'call'

    gk= GarmanKohlhagen(S= spot, K= strike, V= vol, rd= rd, rf= rf, T= T, opt= opt)

    print(f'''
    Reference:
    Spot:   {gk.option.S}    Strike: {gk.option.K}
    Vol:    {gk.option.V}    Tau:{round(gk.option.tau,4)}
    R_d:    {gk.option.rd}   R_f: {gk.option.rf}      

    <Option>  

    Type:  {gk.option.type.upper()}
    Value: {round(gk.option.value, 3)}

    <Greeks>

    Delta:      {round(gk.option.delta,3)}
    Gamma:      {round(gk.option.gamma,3)}
    Vega:       {round(gk.option.vega,3)}
    Theta:      {round(gk.option.theta,3)}
    Rho dmst:   {round(gk.option.rho_d,3)}
    Rho frgn:   {round(gk.option.rho_f,3)}
    ''')


    # import inspect
    # print(inspect.getmembers(gk.option))
