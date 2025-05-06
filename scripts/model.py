import torch
from torch.autograd import Function
from torch import nn
from torch.nn import functional as F

class ReverseLayerF(Function):
    @staticmethod
    def forward(ctx, x, alpha=1):
        ctx.alpha = alpha
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        output = grad_output.neg()*ctx.alpha
        return output, None


class CNNBlock(nn.Module):
    """
    CNN block
    if shortcuts is true, it becomes a resnet block
    """
    def __init__(self,filters,kernel_size=3,shortcuts=True):
        super(CNNBlock,self).__init__()
        padding = int((kernel_size-1)/2)
        self.conv1 = nn.Conv1d(filters, filters, kernel_size=kernel_size, stride=1, padding=padding, bias=False)
        self.bn1   = nn.BatchNorm1d(filters, track_running_stats=True)
        self.conv2 = nn.Conv1d(filters, filters, kernel_size=kernel_size, stride=1, padding=padding, bias=False)
        self.bn2   = nn.BatchNorm1d(filters, track_running_stats=True)
        self.register_buffer('shortcuts',torch.tensor(shortcuts))
    def forward(self,x):
        z = self.conv1(x)        #(,C,L)              
        z = self.bn1(z).relu()
        z = self.conv2(z)        
        z = self.bn2(z)       #(,C,L)            
        # Shortcut connection
        if self.shortcuts:
            z = z + x
        z = z.relu()
        return z
    
class CNNClassifier(nn.Module):
    def __init__(self, n_domains=None, shortcuts=True, n_channels=64, kernel_size = 3, n_blocks = 10):
        super().__init__()
        self.shortcuts = shortcuts
        self.n_channels = n_channels
        self.convIn = nn.Conv1d(4, n_channels, kernel_size=5, stride=1, padding=2, bias=True) 
        self.bnIn   = nn.BatchNorm1d(n_channels, track_running_stats=True)
        self.resblocks = nn.ModuleList([CNNBlock(n_channels, kernel_size, shortcuts=True) for i in range(n_blocks)])
        self.avgpool = nn.AdaptiveAvgPool1d(1)   # pooling each channel of arbitrary size to (1,1)
        self.fcOut   = nn.Linear(n_channels, 2, bias=True)
        if n_domains is not None:
            self.daOut  = nn.Linear(n_channels, n_domains, bias=True)        
        else:
            self.daOut = None
        self.softmax = nn.LogSoftmax(dim=-1)

    def forward(self, x, alpha = 1):
        z = self.convIn(x)
        z = self.bnIn(z)
        z = z.relu()
        for resblock in self.resblocks:
            z = resblock(z)
        z = self.avgpool(z)
        z = z.view(z.size(0), -1)
        # here comes the embedding
        zc = self.fcOut(z)
        if self.daOut is not None:
            zd = ReverseLayerF.apply(z, alpha)
            zd = self.daOut(zd)
            return zc, zd
        else:
            return zc

