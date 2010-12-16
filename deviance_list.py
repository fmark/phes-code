from blist import sortedlist
import operator

class deviance_list:

    def __init__(self):
        self.mean =  0.0
        self._old_mean = 0.0
        self._sum =  0L
        self._n =  0  #n items
        # items greater than the mean
        self._toplist =  sortedlist()
        # items less than the mean
        self._bottomlist = sortedlist(key = operator.neg)
        # Since all items in the "eq list" have the same value (self.mean) we don't need
        # to maintain an eq list, only a count
        self._eqlistlen = 0
        
        self._top_deviance =  0
        self._bottom_deviance =  0

    @property
    def absolute_deviance(self):
        return self._top_deviance + self._bottom_deviance
        
    def append(self,  n):
        # Update summary stats
        self._sum += n
        self._n +=  1
        self._old_mean =  self.mean
        self.mean =  self._sum /  float(self._n)

        # Move existing things around
        going_up = self.mean > self._old_mean
        self._rebalance(going_up)
        
        # Add new item to appropriate list
        if n >  self.mean:
            self._toplist.add(n)
            self._top_deviance +=  n -  self.mean
        elif n == self.mean: 
            self._eqlistlen += 1
        else:
            self._bottomlist.add(n)
            self._bottom_deviance += self.mean -  n
            
            
    def _move_eqs(self,  going_up):
        if going_up:
            self._bottomlist.update([self._old_mean] *  self._eqlistlen)
            self._bottom_deviance += (self.mean - self._old_mean) * self._eqlistlen
            self._eqlistlen = 0
        else:
            self._toplist.update([self._old_mean] *  self._eqlistlen)
            self._top_deviance += (self._old_mean - self.mean) * self._eqlistlen
            self._eqlistlen = 0
        
            
    def _rebalance(self, going_up):
        move_count,  eq_move_count = 0, 0
        if going_up:
            # increase the bottom deviance of the items already in the bottomlist
            if self.mean !=  self._old_mean:
                self._bottom_deviance += len(self._bottomlist) *  (self.mean -  self._old_mean)
                self._move_eqs(going_up)


            # transfer items from top to bottom (or eq) list, and change the deviances
            for n in iter(self._toplist):
                if n < self.mean:
                    self._top_deviance -= n -  self._old_mean
                    self._bottom_deviance += (self.mean -  n)
                    # we increment movecount and move them after the list
                    # has finished iterating so we don't modify the list during iteration
                    move_count +=  1
                elif n == self.mean:
                    self._top_deviance -= n -  self._old_mean
                    self._eqlistlen += 1
                    eq_move_count +=  1
                else:
                    break
            for _ in xrange(0,  move_count):
                self._bottomlist.add(self._toplist.pop(0))
            for _ in xrange(0,  eq_move_count):
                self._toplist.pop(0)

            # decrease the top deviance of the items remain in the toplist
            self._top_deviance -= len(self._toplist) *  (self.mean -  self._old_mean)
        else:
            if self.mean !=  self._old_mean:
                self._top_deviance += len(self._toplist) *  (self._old_mean -  self.mean)
                self._move_eqs(going_up)
            for n in iter(self._bottomlist): 
                if n > self.mean:
                    self._bottom_deviance -= self._old_mean -  n
                    self._top_deviance += n -  self.mean
                    move_count += 1
                elif n == self.mean:
                    self._bottom_deviance -= self._old_mean -  n
                    self._eqlistlen += 1
                    eq_move_count +=  1
                else:
                    break
            for _ in xrange(0,  move_count):
                    self._toplist.add(self._bottomlist.pop(0))
            for _ in xrange(0,  eq_move_count):
                self._bottomlist.pop(0)

            # decrease the bottom deviance of the items remain in the bottomlist
            self._bottom_deviance -= len(self._bottomlist) *  (self._old_mean -  self.mean)

                    
if __name__ ==  "__main__":
    import random
    dv =  deviance_list()
    # Test against some random data,  and calculate result manually to ensure correctness
    rands = [random.randint(-100,  100) for _ in range(0,  100)]
    ns = []
    for n in rands: 
        dv.append(n)
        ns.append(n)
        print("added:%4d,  mean:%3.2f,  oldmean:%3.2f,  mean ad:%3.2f" %
              (n, dv.mean,  dv._old_mean,  dv.absolute_deviance / dv.mean))
        assert sum(ns) == dv._sum,  "Sums not equal!"
        assert len(ns) == dv._n,  "Counts not equal!"
        m = sum(ns) / float(len(ns))
        assert m == dv.mean,  "Means not equal!"
        real_abs_dev = sum([abs(m - x) for x in ns])
        assert abs(real_abs_dev - dv.absolute_deviance) < 0.01, (
            "Absolute deviances not equal. Real:%.2f,  calc:%.2f" %  (real_abs_dev,  dv.absolute_deviance))
    # time it...
    import time
    rands = [random.randint(-100,  100) for _ in range(0,  10000)]
    dv =  deviance_list()
    t =  time.time()
    for n in rands: 
        dv.append(n)
        ad =  dv.absolute_deviance
    tt =  time.time() -  t
    print "deviant list took %.2fs on % d items" %  (tt,  len(rands))
    
    ns =  []
    t =  time.time()
    for n in rands: 
        ns.append(n)
        ad =  sum([abs(m - x) for x in ns])
    tt =  time.time() -  t
    print "standard list took %.2fs on % d items" %  (tt,  len(rands))



    
    
