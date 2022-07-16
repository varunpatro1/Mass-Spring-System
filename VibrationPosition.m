function x = VibrationPosition(x0,m,k,c,f,dt,type)

% stores given information into below variables
x_initial = x0(1);
v_initial = x0(2);


% for ease of calculations below
omega = sqrt(k/m);
xi = (c)/(2*sqrt(m*k));
a = 2.5;

% Forward Ruler
if type == 1  
    x_next = x_initial + (dt*v_initial);
    v_next = v_initial + (dt*((-2*xi*omega*v_initial)-((power(omega,2))*x_initial)+f));
    
    x_initial = x_next;
    v_initial = v_next;
    
    x = [x_initial, v_initial];

% RK-2
elseif type == 2 
    cx1 = dt*v_initial;
    cv1 = dt*((-2*xi*omega*v_initial)-((power(omega,2))*x_initial)+f);
    cx2 = dt*(v_initial+(0.5*cv1));
    cv2 = dt*((-2*xi*omega*(v_initial+(0.5*cv1)))-((power(omega,2))*(x_initial+(0.5*cx1)))+f);

    x_next = x_initial + cx2;
    v_next = v_initial + cv2;
    
    x_initial = x_next;
    v_initial = v_next;
    
    x = [x_initial, v_initial];

% RK-4
elseif type == 4
    cx1 = dt*v_initial;
    cv1 = dt*((-2*xi*omega*v_initial)-((power(omega,2))*x_initial)+f);
    cx2 = dt*(v_initial+(0.5*cv1));
    cv2 = dt*((-2*xi*omega*(v_initial+(0.5*cv1)))-((power(omega,2))*(x_initial+(0.5*cx1)))+f);    
    cx3 = dt*(v_initial+(0.5*cv2));
    cv3 = dt*((-2*xi*omega*(v_initial+(0.5*cv2)))-((power(omega,2))*(x_initial+(0.5*cx2)))+f);
    cx4 = dt*(v_initial+cv3);
    cv4 = dt*((-2*xi*omega*(v_initial+(0.5*cv3)))-((power(omega,2))*(x_initial+(0.5*cx3)))+f);
    
    x_next = x_initial+((1/6)*(cx1+(2*cx2)+(2*cx3)+cx4));
    v_next = v_initial+((1/6)*(cv1+(2*cv2)+(2*cv3)+cv4));
    
    x_initial = x_next;
    v_initial = v_next; 
    
    x = [x_initial, v_initial];
    
else
end