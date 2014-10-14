#! /usr/bin/ruby

class XYZ
  def initialize
    @x = Array.new
    @y = Array.new
    @z = Array.new
    @wl = Array.new
    File.open("xyz.txt").each do |line|
      a = line.chomp.split(/\t/)
      if (a[0].to_f > 0)
        # data
        @wl << a[0].to_f
        @x << a[1].to_f
        @y << a[2].to_f
        @z << a[3].to_f
      end
    end
    @n = @wl.size
  end

  def get(func)
    xx = yy = zz = 0
    (0..(@n-1)).each do |i|
      xx += @x[i]*func.f(@wl[i])
      yy += @y[i]*func.f(@wl[i])
      zz += @z[i]*func.f(@wl[i])
    end
    [xx,yy,zz]
  end
  
  def get_xyz
    xx = yy = zz = 0
    (0..(@n-1)).each do |i|
      xx += @x[i]*yield(@wl[i])
      yy += @y[i]*yield(@wl[i])
      zz += @z[i]*yield(@wl[i])
    end
    [xx,yy,zz]
  end

  def XYZ.xy(xyz)
    v = xyz[0]+xyz[1]+xyz[2]
    [xyz[0]/v, xyz[1]/v]
  end

end

class Planck
  attr_accessor :temp
  @@hckb = 0.0143877696

  def initialize(temp)
    @temp = temp
  end

  def f(wl)
    1.0/wl**5/(Math.exp(@@hckb/wl*1e9/@temp) -1)
  end
end

class XYZ
  def get_black_emitter(temp)
    self.get(Planck.new(temp))
  end
end

class Gaussian
  def initialize(center, std)
    @center = center
    @std = std
  end
  
  def f(wl)
    1.0/Math.sqrt(2*Math::PI)/@std*Math.exp(-(wl-@center)**2/2.0/@std**2)
  end
end
