function [c,t]= poly_coeffs(fcn,var_str)
if nargin==1
    var=symvar(fcn);
    if length(var)>1
       error('Specify var please.');
    end
else
    var=sym(var_str);
end
if isempty(var)
    c0=fcn;
    t0=sym(1);
    disp('Note that: var=[]');
else
    [c0,t0]=coeffs(fcn,var);
end
if (length(t0)==1)&&(t0==1)
    c=c0;
    t.pwr=0;
    t.var=var;
else   
    p_str=strrep(char(t0),'matrix([','');   
    p_str=strrep(p_str,']])',']');         
    p_str=strrep(p_str,'1]','0]');         
    if nargin>=2
         p_str=strrep(p_str,[var_str,'^'],'');  
         p_str=strrep(p_str,var_str,'1');
    else        
         p_str=strrep(p_str,[char(var),'^'],'');  
         p_str=strrep(p_str,char(var),'1');   
    end  
    pwr=eval(p_str);         
     
    m=pwr(1);                              
    c=sym(zeros(1,m+1));                  
    if nargout==2   
        c(m+1-pwr)=c0;
        t.pwr=m:-1:0;
        t.var=var;
    else
        c(pwr+1)=c0;
    end
end