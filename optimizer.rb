
class Symplex_Optimizer_Multi
  attr_accessor :n, :stacked
  def initialize
    @n = 20
    @stacked = 5
  end

  def optimize(param,min=nil,max=nil, &function)
    start_param = param.dup
    best_param = nil
    best_val = nil
    prev_val = 0.0
    count = 0
    @n.times do
      param = Symplex_Optimizer.new.optimize(param, min, max){|p| function.yield(p)}
      val = function.yield(param)
      puts [val, param].flatten.join("\t")
      if (best_param == nil)
        best_param = param
        best_val = val
      end
      if (best_val > val)
        best_param = param
        best_val = val
      end
      if ( (prev_val - val).abs < 1e-6)
        count+=1
        if (count > @stacked)
          param = start_param.map{|x| x+10*(2*rand-1)}
          count = 0
        end
      end
      prev_val = val
    end
    best_param
  end
end

class Symplex_Optimizer
  def initialize
    @limit = 5
    @step = 200
    @tol_cond = 0.0001
    @max_default = 1000
    @min_default = 0
  end

  def prepare(function, param)
    @points = Array.new
    @points << [function.yield(param), param]
    param.size.times do
      new_param = param.zip(@min,@max).map { |x,a_min,a_max| v = x*@limit*rand + (a_max-a_min)*0.8*rand; v = a_min if v < a_min; v = a_max if a_max < v; v }
      @points << [function.yield(new_param), new_param]
    end
  end

  def optimize(param, min=nil, max=nil, &function)
    # function.calc(param)    
    # assumes params are non-negative numbers
    if (min==nil)
      @min = param.map{|x| @min_default}
    else
      @min = min
    end
    if (max==nil)
      @max = param.map{|x| @max_default}
    else
      @max = max
    end
    prepare(function, param)
    n = param.size
    1000.times do
#    while (true)
      @points.sort!{|a,b| a[0]<=>b[0]}
      rtol = 2.0*(@points[0][0] - @points[-1][0])/
        (@points[0][0] + @points[-1][0])
      # if best and worst value are close, return best value
      if (rtol.abs < @tol_cond)
        return @points[0][1]
      end
      
      @points = try(function, @points, -1.0)
#      puts @points[0].flatten.join("\t")
#      puts @points[-1].flatten.join("\t")
      if (@points.last[0] < @points.first[0])
        @points = try(function, @points, 2.0)
      elsif (@points.last[0] >= @points[-2][0])
        y_save = @points.last[0]
        @points = try(function, @points, 0.5)
        if (@points.last[0] >= y_save)
#          puts "shrinking"
#          puts @points[0].flatten.join("\t")
          (1..@points.size-1).each do |i|
#            puts @points[i].flatten.join("\t")
            (0..n-1).each do |j| 
              @points[i][1][j] = 0.5*(@points[0][1][j] + @points[i][1][j])
            end         
            @points[i][0] = function.yield(@points[i][1])
#            puts @points[i].flatten.join("\t")
          end
        end
      end
    end
    return @points[0][1]
  end

  def try(function, points, scale)
    # assume points are sorted
    # points.first : best
    # points.last : worst
    y_worst = points.last[0]
    p_worst = points.last[1]
    n = p_worst.size
    s1 = (1.0-scale)/n
    s2 = s1 - scale

    p_sum = get_sum(points)
    p_try = Array.new
    p_worst.each_index do |i|
      p_try[i] = p_sum[i]*s1-p_worst[i]*s2
      p_try[i] = @min[i] if (p_try[i] < @min[i])
      p_try[i] = @max[i] if (p_try[i] > @max[i])
    end
    y_try = function.yield(p_try)
    if (y_try < y_worst)
      points[-1] = [y_try, p_try]
    end
    return points
  end

  def get_sum(points)
    sum = nil
    points.each do |y|
      point = y[1]
      if (sum == nil)
        sum = point.dup
      else
        sum.each_index do |i|
          sum[i] += point[i]
        end
      end
    end
    sum
  end
end

class MonteCarlo_Optimizer
  def optimize(param, &function)
    # function.calc(param)
    param = param.dup
    n_attempt = 100
    n_temp = 5
    n_high = 100
    n_low = 100

    beta = 1.0/function.yield(param)*0.05
    beta_low = beta*100
    
    tick_p = 0.1
    tick_b = 1.5
    
    result = Array.new
    
    current_val = function.yield(param)
    n_attempt.times do
      n_temp.times do 
        n_high.times do
          trial = param.map{|x| x*(1 + tick_p*(2*rand-1.0))}
          trial_val = function.yield(trial)
          puts [param,trial_val].flatten.join("\t")
          if (Math.exp(-(trial_val - current_val)*beta) > rand)
            # move accepted
            param = trial
          end
        end
        beta *= (1 + tick_b)
        puts beta
        current_val = function.yield(param)
      end
      beta = beta_low
      trial = param.map{|x| x*(1 + tick_p*(2*rand-1.0))}
      trial_val = function.yield(trial)
      puts [param,trial_val].flatten.join("\t")
      if (Math.exp(-(trial_val - current_val)*beta) > rand)
        # move accepted
        param = trial
      end      
      result << [trial_val, param]
    end
    result.sort{|a,b| a[0]<=>b[0]}.first[1]
  end
end

class Golden_Optimizer
  @@tiny = 1e-20
  @@gold = 1.618034
  @@glimit = 100
  @@dx = 1.0
  @@r = 0.61803399
  @@c = 1.0-@@r
  attr_accessor :tol

  def initialize
    @tol = 1e-6
  end

  # function.calc(param)
  def optimize(initial_val, &function)
    #find range
    (a,b,c) = estimate_minima_range(function, initial_val, initial_val+@@dx)
#    puts [a,b,c]

    x0 = x1 = x2 = x3 = 0
    x0 = a
    x3 = c
    if ( (c-b).abs > (b-a).abs)
      x1 = b
      x2 = b+@@c*(c-b)
    else
      x2 = b
      x1 = b-@@c*(b-a)
    end

    (f0, f1, f2, f3) = [x0,x1,x2,x3].map{|x| function.yield(x)}
    while ( (x3-x0).abs > @tol*(x1.abs+x2.abs))
#      puts [x0, x1, x2, x3].join("\t")
      if (f2 < f1)
        (x0, x1, x2) = x1, x2, @@r*x2+@@c*x3
        (f1, f2) = f2, function.yield(x2)
      else
        (x3, x2, x1) = x2, x1, @@r*x1+@@c*x0
        (f2, f1) = f1, function.yield(x1)
      end
    end
    if (f1 < f2)
      return [x1, f1]
    else
      return [x2, f2]
    end
  end

  # code from NR
  def estimate_minima_range(function, a, b)
    fa = function.yield(a)
    fb = function.yield(b)
    if (fb > fa)
      (a,b) = b,a
      (fa, fb) = fb,fa
    end
    c = b+@@gold*(b-a)
    fc = function.yield(c)
    u = fu = 0
    while (fb > fc)
      r = (b-a)*(fb-fc)
      q = (b-c)*(fb-fa)
      denom = (q-r).abs > @@tiny ? q-r : @@tiny
      u = b-( (b-c)*q-(b-a)*r )/ (2.0*denom)

      ulim = b+@@glimit*(c-b)

      if ( (b-u)*(u-c) > 0)
        # u between b, c
        fu = function.yield(u)
        return [b,u,c] if (fu < fc) # min is between a, c
        return [a,b,u] if (fu > fb) # min is between a, u
        u = c + @@gold*(c-b)
        fu = function.yield(u)
      elsif ( (c-u)*(u-ulim) > 0)
        fu = function.yield(u)
        if (fu < fc)
          (b,c,u) = c, u, u+@@gold*(u-c)
          (fb, fc, fu) = fc, fu, function.yield(u)
        end
      elsif ( (u - ulim)*(ulim-c) > 0)
        u = ulim
        fu = function.yield(u)
      else
        u = c+@@gold*(c-b)
        fu = function.yield(u)
      end
      (a,b,c) = b,c,u
      (fa,fb,fc) = fb, fc, fu
    end
    return [a,b,c]
  end
end

class Brent_Optimizer < Golden_Optimizer
  @@iteration_max = 100
  @@zeps = 1e-10
  @@cgold = 0.3819660
  @@dx = 1.0

  def initialize()
    @tol = 1e-6
  end

  def optimize(initial_val, &function)
    #find range
    (a,b,c) = estimate_minima_range(function, initial_val, initial_val+@@dx)
    aa = a < c ? a : c
    bb = a > c ? a : c
    x=w=v=b
    fw=fv=fx=function.yield(b)
    e = tol1 = tol2 = d = 0
    
    @@iteration_max.times do 
      xm = 0.5*(aa+bb)
      tol1 = @tol*x.abs+@@zeps
      tol2 = 2.0*tol1
      if ( (x-xm).abs <= tol2-0.5*(bb-aa))
        return [x, fx]
      end
      if (e.abs > tol1)
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        p = (x-v)*q-(x-w)*r
        q=2.0*(q-r)
        p = -p if q > 0.0
        q = q.abs
        etemp = e
        e = d
        if (p.abs >= (0.5*q*etemp).abs || p <= q*(aa-x) || p >= q*(bb-x))
          e = x>=xm ? aa-x : bb-x
          d = @@cgold*e
        else
          d = p/q
          u = x+d
          d = xm-x>0 ? tol1.abs : -tol1.abs if (u-aa < tol2 || bb-u < tol2)
        end
      else
        e = x>=xm ? aa-x : bb-x
        d = @@cgold*e
      end
      u = d.abs >= tol1 ? x+d : x + (d > 0 ? tol1: -tol1)
      fu = function.yield(u)
      if (fu <= fx)
        if (u >= x)
          aa = x
        else
          bb = x
        end
        (v,w,x) = w,x,u
        (fv,fw,fx) = fw,fx,fu
      else
        if (u < x)
          aa = u
        else
          bb = u
        end
        if (fu <= fw || w == x)
          v = w
          w = u
          fv = fw
          fw = fu
        elsif (fu <= fv || v == x || v == w)
          v = u
          fv = fu
        end
      end
    end
    return [x, fx]
  end
end

class Linearizer
  attr_accessor :pos, :n
  
  def initialize(func, pos, n)
    @func = func
    @pos = pos
    @n = n
  end

  def yield(x)
    param = pos.zip(n).map{|pi,ni| pi+ni*x}
    @func.yield(param)
  end
end

class Powell_Optimizer
  attr_accessor :lin_opt, :lin_opt_tol

  def initialize
    @max_iteration = 200
    @ftol = 1e-4
    @lin_opt = Brent_Optimizer.new
    @lin_opt_tol = 1e-4
  end

  def optimize(param, &func)
    @lin_opt.tol = @lin_opt_tol
    pt = param.dup
    xi = Array.new
    (0..param.size-1).each do |i|
      x = Array.new(param.size, 0)
      x[i]=1
      xi<<x  
    end

    p = param.dup
    fval = func.yield(p)
    lin_func = Linearizer.new(func, p, xi[0])

    @max_iteration.times do 
      fp = fval
      ibig = 0
      del = 0.0
      (0..param.size-1).each do |i|
        fptt = fval
        lin_func.pos = p
        lin_func.n = xi[i]
        (min_x, fval) = @lin_opt.optimize(0){|x| lin_func.yield(x)}
        p = p.zip(xi[i]).map{ |pi,ni| pi+ni*min_x}
        if ( (fptt - fval) > del )
          del = (fptt-fval).abs
          ibig = i
        end
      end
      return [p, fval] if (2.0*(fp - fval).abs <= @ftol*(fp.abs + fval.abs))

      ptt = p.zip(pt).map{|pi,pti| 2.0*pi-pti}
      xit = p.zip(pt).map{|pi,pti| pi-pti}
      pt = p.dup
      fptt = func.yield(ptt)
      if (fptt < fp)
        t = 2.0*(fp-2.0*fval+fptt)*(fp-fval-del)**2 - del*(fp-fptt)**2
        if (t < 0)
          lin_func.pos = p
          lin_func.n = xit
          (min_x, fval) = @lin_opt.optimize(0){|x| lin_func.yield(x)}
          p = p.zip(xit).map{ |pi,ni| pi+ni*min_x}
          xi[ibig] = xi[-1]
          xi[-1] = xit
        end
      end
    end
    return [param,-1]
  end
end

class A_Func3
  attr_reader :count
  def initialize
    @count=0
  end

  def calc(param)
    @count+=1
    param[0]**2+param[1]**2+10*Math.cos(param[0]*2)*Math.cos(param[1]*2)
  end
end

def powell_test()
  opt = Powell_Optimizer.new
  puts "brent"
  opt.lin_opt = Brent_Optimizer.new
  opt.lin_opt_tol = 1e-2
  10.times do
    pos = [rand*10, rand*10]
    f = A_Func3.new
    puts [pos, opt.optimize(pos){|x| f.calc(x)}, f.count].flatten.join("\t")
  end
  puts "golden"
  opt.lin_opt = Golden_Optimizer.new
  opt.lin_opt_tol = 1e-2
  10.times do
    pos = [rand*10, rand*10]
    f = A_Func3.new
    puts [pos, opt.optimize(pos){|x| f.calc(x)}, f.count].flatten.join("\t")
  end
end

class Barrier_Optimizer
  attr_accessor :opt
  def initialize
    @opt = Powell_Optimizer.new
    @opt.lin_opt = Brent_Optimizer.new
  end
  
  def range_penalty(p, min, max)
#    puts [p,min,max, (p-min>0?0:min-p) + (max-p>0?0:max-p)].join("\t")
    (p-min>0?0:min-p) + (max-p>0?0:p-max)
  end

  def optimize(ini_param, min_param, max_param, &func)
    current_param=ini_param
    val=0
    [1, 1e2, 1e4, 1e8].each do |penalty|
      (current_param, val) = @opt.optimize(current_param) do |x|
        func.yield(x) + penalty*x.zip(min_param, max_param).inject(0){|sum, a| sum+range_penalty(a[0], a[1], a[2])}
      end
    end
    [current_param, val]
  end
end

def barrier_test
  puts "barrier_test"
  10.times do
    pos = [rand*10, rand*10]
    puts [pos, Barrier_Optimizer.new.optimize(pos,[1,-1],[1000,1000]){|x| x[0]**2+x[1]**2+10*Math.cos(x[0]*2)*Math.cos(x[1]*2) }].join("\t")
  end
end
#powell_test
#barrier_test

class Mesh_Optimizer
  attr_accessor :opt, :fix, :verbose
  
  def initialize
    @opt = Barrier_Optimizer.new
    @fix = true
    @verbose=false
  end

  def optimize(ini_param, min_param, max_param, mesh, &func)
    best_p = ini_param.dup
    best_val = func.yield(ini_param)
    ini_param_array = (0..(mesh.size-1)).map do |i|
      n = mesh[i]
      if (n>0)
        (1..(n-1)).map{|ii| min_param[i]+(max_param[i]-min_param[i])*ii.to_f/n.to_f}
      else
        [ini_param[i]]
      end
    end
    fix_array = mesh.map{|n| @fix&&n>0?true:false}

    ini_params = ini_param_array.shift.product(*ini_param_array)
    ini_params.each do |test_ini_val|
      tmp_min = min_param.zip(test_ini_val, fix_array).map{|a| a[2]?a[1]:a[0]}
      tmp_max = max_param.zip(test_ini_val, fix_array).map{|a| a[2]?a[1]:a[0]}
      (tmp_p, tmp_val) = @opt.optimize(test_ini_val, tmp_min, tmp_max){|x| func.yield(x)}
      puts [test_ini_val, tmp_p, tmp_val].flatten.join("\t") if (verbose)
      if (tmp_val < best_val)
        best_val = tmp_val
        best_p = tmp_p
      end
    end
#    puts [best_p, best_val].join("\t")
    @opt.optimize(best_p, min_param, max_param){|x| func.yield(x)}
  end
end

#puts [Mesh_Optimizer.new.optimize([0,0],[0,0],[10,10],[10,0]){|x| x[0]**2+x[1]**2+10*Math.cos(x[0]*2)*Math.cos(x[1]*2) }].join("\t")

