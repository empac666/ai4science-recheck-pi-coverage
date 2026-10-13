import unittest
from dataclasses import replace
import numpy as np
from nextrun.core import (feature, euler_step, trajectory, observe, indices, noise_scales, frozen,
                          fit, Model, predict, history, initial_history, reveal, quantile_width,
                          halfwidth, score, evaluate, A_METHODS)
from nextrun.statistics import (bias_diagnostic, describe, paired, ratio, compatibility,
                                raw_m, raw_e, raw_j, gate, conclusion, Q)


class ProtocolChecks(unittest.TestCase):
    def test_features_order_and_penalty(self):
        x=np.arange(1.,21.)
        expected=[1.]+list(x)+[x[a]*x[b] for a in range(20) for b in range(a,20)]
        np.testing.assert_array_equal(feature(x),expected)
        self.assertEqual(len(expected),231)
        # Nonformal, artificial fixed design; no preregistered seed is used.
        rng=np.random.Generator(np.random.PCG64(17))
        lags=rng.normal(size=(1600,20));y=0.3+lags[:,0]*2+rng.normal(size=1600)
        model=fit(lags,y);r=feature(lags)
        np.testing.assert_allclose((r.T@r+16*np.eye(231))@model.w,r.T@y,atol=1e-10)
        self.assertAlmostEqual(model.tau,np.sqrt(np.mean((y-r@model.w)**2)))
        self.assertFalse(model.w.flags.writeable)
        self.assertTrue(0<=model.df<=231)
        self.assertTrue(1<=model.tau_df/model.tau<=40/37)
        # Cal/test data are absent from the fitting API and cannot enter its Gram/RHS.
        calibration=np.full((800,20),999.)
        before=model.w.copy();feature(calibration)@model.w
        np.testing.assert_array_equal(model.w,before)

    def test_synchronous_euler_and_first_sample(self):
        x=np.array([2.,5.,7.]);t=.3;h=.005
        expected=x+h*np.array([np.sin(t)-.1*2+.5*np.tanh(5),2-.2*5-2*np.pi*7,2*np.pi*5-.2*7])
        np.testing.assert_array_equal(euler_step(x,t,h),expected)
        state=np.zeros(3)
        for k in range(20): state=euler_step(state,k*.005,.005)
        np.testing.assert_array_equal(trajectory(.005,burn=0,count=1)[0],state)
        np.testing.assert_array_equal(trajectory(.005,burn=1,count=1)[0],trajectory(.005,burn=0,count=2)[1])

    def test_observation_and_shift_indices(self):
        for s,want in [(0,2),(.15,2.3),(.35,2.7)]:
            self.assertEqual(observe(np.array([2,5,7]),2,s),want)
            self.assertEqual(observe(np.array([2,50,70]),2,s),want)
        for off in (20,21):
            inputs,targets=indices(off);sigma=noise_scales(off,'shift')
            self.assertEqual(targets[0],off);self.assertEqual(targets[-1],3199+off)
            self.assertEqual(sigma[2400+off+159],.15);self.assertEqual(sigma[2400+off+160],.35)
            for t,n in [(160,0),(179,19),(180,20)]:
                self.assertEqual(np.count_nonzero(sigma[inputs[2400+t]]==.35),n)
            self.assertTrue(np.all(inputs<targets[:,None]))
        i,q=indices();self.assertEqual(q[2560],2581);self.assertEqual(i[0,-1],1)

    def test_snapshot_isolation_and_history_ids(self):
        base=np.arange(-800.,0.);state=initial_history(base,'A');old=state.values.copy()
        with self.assertRaises(ValueError): state.values[0]=0
        base[-1]=999
        np.testing.assert_array_equal(state.values,old)
        state2=reveal(state,'A',-3,True,0)
        np.testing.assert_array_equal(state.values,old)
        self.assertEqual(state2.values[-1],3)
        self.assertFalse(np.shares_memory(state.values,state2.values))
        other=initial_history(base,'A');self.assertFalse(np.shares_memory(other.values,state.values))
        for method in ('G','A'):
            s=initial_history(np.arange(800.),method)
            for t in range(261):
                if t in (160,180,260): np.testing.assert_array_equal(s.ids,np.arange(t-100,t))
                s=reveal(s,method,(-1)**t,False,t)
        full=initial_history(np.arange(-800.,0.),'A_full')
        self.assertEqual(len(reveal(full,'A_full',-1,False,0).values),801)

    def test_adaptive_boundaries_and_order_statistics(self):
        h=history(np.arange(1.,101.),np.arange(-100,0),absolute=True)
        self.assertAlmostEqual(reveal(h,'A',-1,True,0).alpha,.0405)
        self.assertAlmostEqual(reveal(h,'A',-1,False,0).alpha,.0505)
        self.assertGreaterEqual(quantile_width(h.values,.04,'A'),quantile_width(h.values,.05,'A'))
        self.assertEqual(quantile_width(h.values,0,'A'),100)
        self.assertEqual(quantile_width(h.values,1,'A'),1)
        self.assertEqual(quantile_width(h.values,.05,'A_order'),96)
        self.assertEqual(quantile_width(h.values,0,'A_order'),np.inf)
        self.assertEqual(quantile_width(h.values,1,'A_order'),0)
        self.assertEqual(quantile_width(h.values,0,'A_inf'),np.inf)
        self.assertEqual(quantile_width(h.values,.05,'A_inf'),quantile_width(h.values,.05,'A'))
        for bad in (-.01,1.01,np.nan):
            with self.assertRaises(ValueError): quantile_width(h.values,bad,'A')
        with self.assertRaises(AssertionError): quantile_width(np.array([-1.,2.]),.05,'A')
        zero=replace(h,alpha=0.)
        self.assertEqual(reveal(zero,'A',5,True,0).alpha,0)
        self.assertEqual(reveal(zero,'A_clip',5,True,0).alpha,.001)
        self.assertEqual(np.partition(np.arange(1.,801.),760)[760],761)
        self.assertEqual(halfwidth('rolling_F',None,0,0,h),96)

    def test_interval_score_and_invalid_width(self):
        for y in (1,3): self.assertEqual(score(y,1,3),(True,2,2))
        self.assertEqual(score(0,1,3),(False,2,42))
        self.assertEqual(score(4,1,3),(False,2,42))
        self.assertEqual(score(100,-np.inf,np.inf),(True,np.inf,np.inf))
        for h in (-1,np.inf,np.nan):
            with self.assertRaises(ValueError): halfwidth('F',None,0,h)
        with self.assertRaises(ValueError): score(0,np.nan,1)

    def test_current_future_perturbation_cannot_change_current_interval(self):
        model=Model(frozen(np.zeros(231)),frozen(4*np.eye(231)),.2,1.,.21)
        y=np.sin(np.arange(3221)*.3)*.1
        target=np.arange(2421,3221);cal=np.linspace(-.3,.3,800)
        methods=['B','F','G','A','A_full','A_clip','A_order','A_inf','G_L','rolling_F','B_df']
        a=evaluate(model,y,target,cal,methods)
        changed=y.copy();changed[target[160]:]+=1
        b=evaluate(model,changed,target,cal,methods)
        for key in ('lower','upper','alpha','halfwidth'):
            np.testing.assert_array_equal(a[key][:,:161],b[key][:,:161])
        np.testing.assert_array_equal(a['center'][:161],b['center'][:161])
        self.assertTrue(np.all(a['hit'][8]>=a['hit'][2]))
        frozen_lags=frozen(y[:20]);pred=predict(model,frozen_lags)
        y[:20]=20;self.assertEqual(pred,predict(model,frozen_lags))

    def test_quality_gate_compatibility_and_strict_boundaries(self):
        self.assertEqual(compatibility(.78,[.77,.79],.726,.05),'抽样不确定下仍相容，需披露')
        self.assertEqual(compatibility(.78,[.777,.79],.726,.05),'不相容')
        self.assertEqual(compatibility(.67,[.65,.675],.726,.05),'不相容')
        self.assertNotEqual(compatibility(.78,[.776,.8],.726,.05),'不相容')
        self.assertEqual(compatibility(0,None,.941,.02),'不可估计')
        self.assertEqual(compatibility(.91,[.90,.92],.941,.02),'不相容')
        for lo in (.01,.15): self.assertFalse(raw_m([lo,.2]))
        self.assertTrue(raw_m([.15001,.2]))
        for ci in ([-.02,.01],[-.01,.02]): self.assertFalse(raw_e(ci))
        ids=np.array([[0,0],[0,1],[1,1]])
        d=describe([.7,.8],ids)
        np.testing.assert_array_equal(d['compatibility_ci'],np.quantile([.7,.75,.8],[.05/12,1-.05/12],method='linear'))

    def test_cost_dispatch_bias_and_joint_gate(self):
        good={'estimable':True,'ci90':[-2,-1]};edge={'estimable':True,'ci90':[-1,0]}
        self.assertTrue(raw_j(good,good));self.assertFalse(raw_j(edge,good));self.assertFalse(raw_j(good,edge))
        ids=np.array([[0,1],[1,0]])
        self.assertFalse(paired(np.array([1.,np.inf]),np.ones(2),ids)['estimable'])
        self.assertFalse(describe([1,np.inf],ids)['estimable'])
        self.assertFalse(ratio([1,2],[0,2],ids)['estimable'])
        self.assertFalse(ratio([1,2],[np.inf,2],ids)['estimable'])
        self.assertEqual(bias_diagnostic(1,np.ones(3))['b'],'N/A (degenerate)')
        self.assertEqual(bias_diagnostic(0,np.ones(3))['b'],'ERROR (degenerate non-zero)')
        self.assertTrue(bias_diagnostic(0,np.arange(3.))['bias_flag'])
        # Auxiliary ratio flags cannot change linear J.
        self.assertTrue(raw_j(good,good))
        self.assertEqual(gate(False,True,True,True)['M'],'未检验（质量门关闭）')
        self.assertEqual(gate(True,False,True,True)['J'],'未检验（前序未通过）')
        self.assertEqual(gate(True,True,False,True)['J'],'未检验（前序未通过）')
        self.assertEqual(gate(True,True,True,False,False)['J'],'不可估计')

    def test_conclusion_priority_and_limits(self):
        s=gate(True,True,True,True)
        text=conclusion(s,True,True,True,False)
        self.assertIn('对积分/索引约定敏感',text['M'])
        self.assertNotIn('支持平均覆盖改善',text['M'])
        for piece in ('依赖本稿预先固定的 A','对积分/索引约定敏感','检查未完成'):
            self.assertIn(piece,text['E'])
        self.assertTrue(text['M'].startswith(Q));self.assertTrue(text['E'].startswith(Q))
        self.assertEqual(text['E'].count(Q),1)
        blocked=conclusion(gate(False,True,True,True),False,False,False,True)
        self.assertIsNone(blocked['E']);self.assertIn('主推断被阻断',blocked['M'])
        # N deliberately cannot override the M sensitivity branch.
        self.assertIn('敏感',conclusion(s,True,False,False,True)['M'])

if __name__=='__main__': unittest.main(verbosity=2)

class ReportingChecks(unittest.TestCase):
    def test_early_width_and_bias_in_same_paragraph(self):
        from nextrun.reporting import main_paragraph
        stub={'conclusion':{'M':Q+'主设定下预注册门槛通过，但该结论对积分/索引约定敏感。'},
              'comparisons':{'G-B/shift/early50':{'coverage':{'mean':.1,'se':.01,'ci95':[.08,.12]}}},
              'width_reference':{'B':{'over_10_percent':True,'relative_difference':-.2},'G':{'over_10_percent':False}},
              'J_statistics':{'width':{'D':{'diagnostic':{'bias_flag':False,'b':0.}},'ratio':{'diagnostic':{'bias_flag':True,'b':.2}}}},
              'gated':{'J':'通过'}}
        text=main_paragraph(stub)
        for phrase in ('前 50 点','SE=','95% 配对区间','无确认性资格','覆盖结论不附带宽度复现','bootstrap 偏差相对 SE 超过预定阈值','不参与线性 J'):
            self.assertIn(phrase,text)
        self.assertNotIn('\n',text)
