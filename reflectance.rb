#! /usr/bin/ruby
require 'matrix'
require 'complex'
require 'optimizer'
require 'xyz'

class Param
  attr_accessor :wl, :theta, :pol
  def angle=(angle)
    @theta = angle.to_f/180.0*Math::PI
  end

  def angle
    @theta/Math::PI*180.0
  end
end

class Index
  def initialize(n=1.0)
    @n = n
  end

  def n(p)
    # p : an object contains parameters
    @n
  end
end

class RealIndex
  def initialize(path)
    @n = Array.new
    @wl = Array.new
    read_file(path)
  end

  def read_file(path)
    File.open(path) do |file|
      file.each do |line|
        if (line =~ /^\#/)
          #comment
        elsif ( (a = line.split(/\t/)).size == 3)
          @wl << a[0].to_f
          @n << Complex(a[1].to_f,a[2].to_f)
        end
      end
    end
  end

  def interpolate(p)
    if (p.wl <= @wl[0])
      @n[0]
    elsif (p.wl >= @wl[-1])
      @n[-1]
    else
      left = 0
      right = @wl.size
      while (left + 1 < right)
        pivot = (left+right)/2
        if (@wl[pivot] < p.wl)
          left = pivot
        else
          right = pivot
        end
      end

      x_l = @wl[left]
      y_l = @n[left]
      x_r = @wl[right]
      y_r = @n[right]
      y_l + (p.wl - x_l)/(x_r-x_l)*(y_r - y_l)
    end
  end

  def n(p)
    interpolate(p)
  end
end

class Layer
  def initialize(refraction_index,thickness)
    @n = refraction_index
    @d = thickness
  end

  def n(p)
    @n.n(p).conjugate
  end

  def m(p,n_0)
    wl = p.wl
    theta = p.theta
    cos_theta = cos_theta(p,n_0)
    delta = 2*Math::PI/wl*n(p)*@d*cos_theta
    eta = eta(p,n_0)
    Matrix[[Math.cos(delta), Complex(0,1)*Math.sin(delta)/eta],
            [Complex(0,1)*Math.sin(delta)*eta, Math.cos(delta)]]
  end
  
  def eta(p,n_0)
    pol = p.pol
    if (pol == :s)
      eta = n(p)*cos_theta(p,n_0)
    else
      eta = n(p)/cos_theta(p,n_0)
    end
    eta
  end

  def cos_theta(p,n_0)
    theta = p.theta
    Math.sqrt(1-(n_0/n(p))**2*Math.sin(theta)**2)
  end
end

class Layers
  def initialize(l)
#    @layers = Array.new
    @layers = l
  end
  
  def ref(p)
    if (@layers.size >= 2)    
      first_layer = @layers.first
      def first_layer.n(p)
        @n.n(p).real
      end
      last_layer = @layers.last
      def last_layer.n(p)
        @n.n(p).real
      end
      n_0 = @layers.first.n(p)
      eta_0 = @layers.first.eta(p,n_0)
      eta_m = @layers.last.eta(p,n_0)

      m = Matrix[[1,0],[0,1]]
      (1..@layers.size-2).each do |i|
        layer = @layers[i]
        m = m*layer.m(p,n_0)
      end
      mm = m*Matrix[[1],[eta_m]]
      b = mm[0,0]
      c = mm[1,0]
      rho = (eta_0*b-c)/(eta_0*b+c)
      tau = 2*eta_0/(eta_0*b+c)
      rr = rho.abs**2
      tt = eta_m/eta_0*tau.abs**2
      def first_layer.n(p)
        @n.n(p).conjugate
      end
      def last_layer.n(p)
        @n.n(p).conjugate
      end
      [rr,tt]
    else
      0
    end
  end
end

class Array
  def sum
    sum=0
    self.each do |val|
      sum+=val
    end
    sum
  end
  
  def avg
    sum/self.size
  end
end

class RefUtil
  def RefUtil.angle_scan(layers, parameter, hash=nil)
    old_theta = parameter.theta
    angles = Array.new
    r_s = Array.new
    r_p = Array.new
    start_angle = 0
    end_angle = 90.0
    num = 90
    print = true
    if (hash != nil)
      start_angle = hash[:start_angle] if (hash[:start_angle] != nil)
      end_angle = hash[:end_angle] if (hash[:end_angle] != nil)
      num = hash[:num] if (hash[:num] != nil)
      print = hash[:print?] if (hash[:print?] != nil)
    end
    (0..(num-1)).each do |i|
      angle = start_angle + (end_angle - start_angle)/num*i
      angles << angle
      parameter.angle = angle
      parameter.pol = :s
      r_s << layers.ref(parameter)[0]
      parameter.pol = :p
      r_p << layers.ref(parameter)[0]
    end
    if (print)      
      puts ["angle", "r_s", "r_p"].join("\t")
      angles.each_index {|i| puts [angles[i],r_s[i],r_p[i]].join("\t")}
    end
    parameter.theta = old_theta
    [angles, r_s, r_p]
  end

  def RefUtil.wavelength_scan(layers,parameter,hash=nil)
    wls = Array.new
    r_s = Array.new
    r_p = Array.new
    t_s = Array.new
    t_p = Array.new
    
    start_wl = 400
    end_wl = 800
    num = 400
    print = true
    tag = nil
    if (hash != nil)
      start_wl = hash[:start_wl] if (hash[:start_wl] != nil)
      end_wl = hash[:end_wl] if (hash[:end_wl] != nil)
      num = hash[:num] if (hash[:num] != nil)
      num = (end_wl - start_wl).to_f/hash[:tick] if (hash[:tick] != nil)
      print = hash[:print] if (hash[:print] != nil)
      tag = hash[:tag]
    end
    (0..num).each do |i|
      wl = start_wl + (end_wl - start_wl)/num.to_f*i
      wls << wl
      parameter.wl = wl
      parameter.pol = :s
      r_s << layers.ref(parameter)[0]
      t_s << layers.ref(parameter)[1]
      parameter.pol = :p
      r_p << layers.ref(parameter)[0]
      t_p << layers.ref(parameter)[1]
    end
    if (print)
      if (tag != nil)
        puts ["wavelength", "r_s", "t_s", "r_p", "t_p"].map{|x| x+"_"+tag}.join("\t")
      else
        puts ["wavelength", "r_s", "t_s", "r_p", "t_p"].join("\t")
      end
      wls.each_index {|i| puts [wls[i],r_s[i],t_s[i],r_p[i],t_p[i]].join("\t")}
    end
    [wls, r_s, r_p, t_s, t_p]
  end    
  
  def RefUtil.avg_ref(layers, parameter, hash=nil)
    if (hash == nil)
      hash_child = Hash.new
    else
      hash_child = hash.dup
    end
    hash_child[:print] = false
    (wls, r_s, r_p) = wavelength_scan(layers,parameter,hash_child)
    a = [r_s.avg, r_p.avg]
    if (hash==nil || (hash!=nil && hash[:print]))
      puts "r_s\tr_p"
      puts a.join("\t")
    end
    a
  end

  @@xyz = XYZ.new 
  @@d65 = Planck.new(6500)
  def RefUtil.ref_xyz(layers, parameter, hash=nil)
    if (hash == nil)
      hash_child = Hash.new
    else
      hash_child = hash.dup
    end
    hash_child[:print] = false
    @@xyz.get_xyz do |wl|
      parameter.wl = wl
      parameter.pol = :s
      r_s = layers.ref(parameter)[0]
#      t_s << layers.ref(parameter)[1]
      parameter.pol = :p
      r_p = layers.ref(parameter)[0]
#      t_p << layers.ref(parameter)[1]
      (r_s+r_p)*0.5*@@d65.f(wl)
    end
  end
end

