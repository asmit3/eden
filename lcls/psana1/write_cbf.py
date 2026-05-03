from dxtbx.format.cbf_writer import FullCBFWriter
writer = FullCBFWriter("locator.loc")
for i in range(5):
  writer.write_cbf("example_%d.cbf"%i, index=i)
