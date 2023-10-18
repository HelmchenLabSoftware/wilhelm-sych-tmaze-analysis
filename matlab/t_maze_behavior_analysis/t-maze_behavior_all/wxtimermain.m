%% wxtimermain
 % Main program for the automation of the T-maze
 % JL Alatorre Warren

global tmr
global aa
global bb

n = 1;
aa = 0;
bb = 0;

while n < 1000
  
  tmr = timer('ExecutionMode', 'FixedRate', ...
            'TasksToExecute', 1, ...
            'Period', 0.300, ...
            'TimerFcn', {@wxtimercallback, 1, frameMonochrome, roi01, roi02});

  start(tmr);
  
  if aa == 1
    wxservo('start','open')
    wxservo('start','close')
  end
  
  if bb == 1
    wxservo('right','open')
    wxservo('right','close')
  end
  
  frameMonochrome = LucamCaptureMonochromeFrame(cameraNumber);
  n=n+1;
  display(n)
  
end