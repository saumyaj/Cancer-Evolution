function ret = f(pr,ps,cr,cs,k)

global lambda;
ret = lambda*log(k/(pr+ps+cr+cs));