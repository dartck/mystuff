package TutorialExamples

import Chisel._
import scala.collection.mutable.ArrayBuffer

object C {
  // start execution at this address
  val STARTADDR = 0x00400000

  // memory and datapath sizes
  val WORDLENGTH = 32	// datapath size
  val IADDRZ = 10	// number of bits for instruction memory 
  val DADDRZ = 10	// number of bits for data memory 

  // instruction opcodes
  val OP_RTYPE = UInt(0x00)
  val OP_JR = UInt(0x00)
  val OP_J = UInt(0x02)
  val OP_JAL = UInt(0x03)
  val OP_BEQ = UInt(0x04)
  val OP_BNE = UInt(0x05)
  val OP_ADDI = UInt(0x08)
  val OP_ADDIU = UInt(0x09)
  val OP_ORI = UInt(0x0d)
  val OP_LUI = UInt(0x0f)
  val OP_LW = UInt(0x23)
  val OP_SW = UInt(0x2b)

  // function codes for R-type instructions
  val FUNC_ADDU = UInt(0x21)
  val FUNC_SUBU = UInt(0x23)
  val FUNC_SLTU = UInt(0x2b)
  val FUNC_OR = UInt(0x25)
  val FUNC_JR = UInt(0x08)

  // ALU control
  val ALUCTLZ = 4		// number of ALU control bits
  val ALU_ADD = UInt(0x02)
  val ALU_SUB = UInt(0x06)
  val ALU_SLT = UInt(0x07)
  val ALU_OR = UInt(0x01)
  val ALU_AND = UInt(0x00)
  val ALU_NOR = UInt(0x0C)
  val ALU_LUIB = UInt(0x0D)	// output is Cat(B(15, 0), UInt(0, width=16))

  // Registers
  val REG_LR = 31

}

class Dpath extends Module {
  val io = new Bundle {
    val isWr   = Bool(INPUT)
    val wrAddr = UInt(INPUT, C.IADDRZ)
    val wrData = Bits(INPUT, C.WORDLENGTH)
    val boot   = Bool(INPUT)
    val out    = Bits(OUTPUT, C.WORDLENGTH)
    val pc    = UInt(OUTPUT, C.WORDLENGTH)
    val inst    = UInt(OUTPUT, C.WORDLENGTH)
    val alu_src = Bool(INPUT)
    val alu_op = UInt(INPUT, C.ALUCTLZ)
    val reg_write = Bool(INPUT)
    val reg_dst = Bool(INPUT)
    val mem_toreg = Bool(INPUT)
    val mem_write = Bool(INPUT)
    val jump = Bool(INPUT)
    val jrmux = Bool(INPUT)
    val branch = Bool(INPUT)
    val link_write = Bool(INPUT)
    val zero = Bool(OUTPUT)
    val alu1    = UInt(OUTPUT, C.WORDLENGTH)
    val alu2    = UInt(OUTPUT, C.WORDLENGTH)
    val alures    = UInt(OUTPUT, C.WORDLENGTH)
  }
  val regfile = Mem(Bits(width = C.WORDLENGTH), C.WORDLENGTH)
  val imem = Mem(Bits(width = C.WORDLENGTH), 1 << C.IADDRZ)
  val dmem = Mem(Bits(width = C.WORDLENGTH), 1 << C.DADDRZ)
  val pc   = Reg(init=UInt(C.STARTADDR, C.WORDLENGTH))

  // NOP
  val nop = UInt(0, width = C.WORDLENGTH)

  // create reg blocks for pipeline storage
  // IF/ID
  val ifid_pc = Reg(init = UInt(C.STARTADDR, width = C.WORDLENGTH))	//pc
  val ifid_inst = Reg(init = UInt(0, width = C.WORDLENGTH))		//inst

  // ID/EX
  val idex_pc = Reg(init = UInt(C.STARTADDR, width = C.WORDLENGTH)) 	//pc
  val idex_rs = Reg(init = UInt(0, width = C.WORDLENGTH))		//rs
  val idex_rt = Reg(init = UInt(0, width = C.WORDLENGTH))		//rt
  val idex_op = Reg(init = UInt(0, width = 5))				//op
  val idex_rsi = Reg(init = UInt(0, width = 5))				//rsi
  val idex_rti = Reg(init = UInt(0, width = 5))				//rti
  val idex_rdi = Reg(init = UInt(0, width = 5))				//rdi
  val idex_dstreg = Reg(init = UInt(0, width = 5))			//dstreg
  val idex_sextimm = Reg(init = UInt(0, width = C.WORDLENGTH))		//sextimm
  val idex_jumpad = Reg(init = UInt(0, width = C.WORDLENGTH))		//absjumpa
  val idex_branchad = Reg(init = UInt(0, width = C.WORDLENGTH))		//relbrancha
  val idexsig_alusrc = Reg(init = Bool(false))				//alu_src
  val idexsig_aluop = Reg(init = UInt(width = C.ALUCTLZ))		//alu_op
  val idexsig_memwrite = Reg(init = Bool(false))			//mem_write
  val idexsig_branch = Reg(init = Bool(false))				//branch
  val idexsig_jump = Reg(init = Bool(false))				//jump
  val idexsig_jrmux = Reg(init = Bool(false))				//jrmux
  val idexsig_memtoreg = Reg(init = Bool(false))			//mem_toreg
  val idexsig_regwrite = Reg(init = Bool(false))			//reg_write
  val idexsig_linkwrite = Reg(init = Bool(false))			//link_write

  // EX/MEM
  val exmem_pc = Reg(init = UInt(C.STARTADDR, width = C.WORDLENGTH))	//pc
  val exmem_aluout = Reg(init = UInt(0, width = C.WORDLENGTH))		//alu_out
  val exmem_writedata = Reg(init = UInt(0, width = C.WORDLENGTH))	//writedata
  val exmem_dstreg = Reg(init = UInt(0, width = 5))			//dstreg
  val exmemsig_memwrite = Reg(init = Bool(false))			//mem_write
  val exmemsig_memtoreg = Reg(init = Bool(false))			//mem_toreg
  val exmemsig_regwrite = Reg(init = Bool(false))			//reg_write
  val exmemsig_linkwrite = Reg(init = Bool(false))			//link_write

  // MEM/WB              
  val memwb_pc = Reg(init = UInt(C.STARTADDR, width = C.WORDLENGTH))	//pc        
  val memwb_aluout = Reg(init = UInt(0, width = C.WORDLENGTH))		//alu_out
  val memwb_dmemout = Reg(init = UInt(0, width = C.WORDLENGTH))		//dmem_out
  val memwb_dstreg = Reg(init = UInt(0, width = 5))			//dstreg
  val memwbsig_memtoreg = Reg(init = Bool(false))			//mem_toreg
  val memwbsig_regwrite = Reg(init = Bool(false))			//reg_write
  val memwbsig_linkwrite = Reg(init = Bool(false))			//link_write

  // INSTRUCTION FETCH 
  // next stage
  val instruction = imem(pc(C.IADDRZ-1, 2))

  // INSTRUCTION DECODE 
  // decode
  val inst = ifid_inst
  val op   = ifid_inst(31,26)
  val rsi  = ifid_inst(25,21)
  val rti  = ifid_inst(20,16)
  val rdi  = ifid_inst(15,11)
  val shamt  = ifid_inst(10,6)
  val func  = ifid_inst(5,0)
  val imm  = ifid_inst(15,0)
  //val rs = Mux(rsi === Bits(0), Bits(0), regfile(rsi))	// regfile(0) = 0
  //val rt = Mux(rti === Bits(0), Bits(0), regfile(rti))

  // immediate values extracted from instruction
  val sextimm = Cat(Fill(16, imm(15)), imm)
  val absjumpa = Cat(ifid_pc(31, 28), inst(25, 0), UInt(0, width=2))
  val relbrancha = Cat(ifid_pc(C.WORDLENGTH-1, 2) + sextimm - UInt(1), UInt(0, width=2))

  // destination register
  val dest = Mux(io.reg_dst, rdi, rti)
  val stall = Bool()

  // store for next stage
  //idex_rs := rs
  //idex_rt := rt
  idex_op := op
  idex_rsi := rsi
  idex_rti := rti
  idex_rdi := rdi
  idex_dstreg := dest
  idex_sextimm := sextimm
  idex_jumpad := absjumpa
  idex_branchad := relbrancha
  idex_pc := ifid_pc
  idexsig_alusrc := io.alu_src
  idexsig_aluop := io.alu_op
  idexsig_regwrite := io.reg_write
  idexsig_memtoreg := io.mem_toreg
  idexsig_memwrite := io.mem_write
  idexsig_branch := io.branch
  idexsig_jump := io.jump
  idexsig_jrmux := io.jrmux
  idexsig_linkwrite := io.link_write
  
  when ((instruction(20,16) === dest || instruction(25,21) === dest) && dest != UInt(0) && io.reg_write){
     stall := Bool(true)
     ifid_inst := nop
  } .elsewhen ((instruction(20,16) === idex_dstreg || instruction(25,21) === idex_dstreg) && idex_dstreg != UInt(0) && idexsig_regwrite){
     stall := Bool(true)
     ifid_inst := nop
  } .elsewhen (io.branch || idexsig_branch || io.jump || idexsig_jump){
     stall := Bool(true)
     ifid_inst := nop
  } .otherwise {
     stall := Bool(false)
     ifid_inst := instruction
     ifid_pc := pc + UInt(4)
  }

  // EXECUTION
  // alu

  val rs = regfile(idex_rsi)
  val rt = regfile(idex_rti)
  val alu_in1 = rs
  val alu_in2 = Mux(idexsig_alusrc, rt, idex_sextimm)
  val alu_out = MuxCase(UInt(0), Array(
              (idexsig_aluop === C.ALU_ADD) -> (alu_in1 + alu_in2).toUInt,
              (idexsig_aluop === C.ALU_SUB) -> (alu_in1 - alu_in2).toUInt,
              (idexsig_aluop === C.ALU_LUIB) -> (Cat(alu_in2(15, 0), UInt(0, width=16))).toUInt,
              (idexsig_aluop === C.ALU_SLT) -> (Mux(rs < rt, UInt(1), UInt(0))).toUInt,
              (idexsig_aluop === C.ALU_OR) -> (alu_in1 | alu_in2).toUInt))
  val zero = (alu_in1 === alu_in2)
  val dobranch = idexsig_branch && ((zero && idex_op === C.OP_BEQ) ||
                  (!zero && idex_op === C.OP_BNE))
  
  // pc next address
  val stallpc = Mux(stall, pc, (pc + UInt(4)))
  val jumpaddr = Mux(idexsig_jrmux, rs, idex_jumpad)
  val braddr = Mux(dobranch, idex_branchad, stallpc)
  val nextpc = Mux(idexsig_jump, jumpaddr, braddr)
 
  // store for next stage
  exmem_pc := idex_pc
  exmem_aluout := alu_out
  exmem_writedata := rt
  exmem_dstreg := idex_dstreg
  exmemsig_memwrite := idexsig_memwrite
  exmemsig_memtoreg := idexsig_memtoreg
  exmemsig_regwrite := idexsig_regwrite 
  exmemsig_linkwrite := idexsig_linkwrite


  // MEMORY
  // data memory
  val dmem_out = dmem(exmem_aluout)

  // store for next stage
  memwb_pc := exmem_pc
  memwb_aluout := exmem_aluout
  memwb_dmemout := dmem_out
  memwb_dstreg := exmem_dstreg 
  memwbsig_memtoreg := exmemsig_memtoreg
  memwbsig_regwrite := exmemsig_regwrite
  memwbsig_linkwrite := exmemsig_linkwrite

  // WRITE BACK
  // changes of state
  when (io.isWr) {
    imem(io.wrAddr) := io.wrData
  } .elsewhen (io.boot) {
    pc := UInt(C.STARTADDR)
  } .otherwise {
    // debug messages
  //   printf("** pc=%x inst=%x op=%x func=%x alu_in1=%x alu_in2=%x dmem_out=%x alu_out=%x braddr=%x, nextpc=%x\n", pc, inst, op, func, alu_in1, alu_in2, dmem_out, alu_out, braddr, nextpc)
    when (memwbsig_regwrite) {
      val writedata = Mux(memwbsig_memtoreg, memwb_dmemout, memwb_aluout)
      val writedatalink = Mux(memwbsig_linkwrite, memwb_pc, writedata)
      val dstreglink = Mux(memwbsig_linkwrite, UInt(C.REG_LR), memwb_dstreg)
      regfile(dstreglink) := writedatalink
      // printf("RF(%d)=%x ALU(%x,%x,%x)\n", dstreglink, writedatalink, alu_in1, alu_op, alu_in2)
    }


    // MEMORY
    when (exmemsig_memwrite) {
      dmem(exmem_aluout) := exmem_writedata
      // printf("DM(%x)=%x\n", alu_out, rt)
    }
    pc := nextpc
  }


  // io for the processor
  io.out := regfile(UInt(25))
  io.pc := pc
  io.inst := inst
  io.zero := zero
  io.alu1 := alu_in1
  io.alu2 := alu_in2
  io.alures := alu_out
}  

class Single extends Module {
  val io = new Bundle {
    val isWr   = Bool(INPUT)
    val wrAddr = UInt(INPUT, C.IADDRZ)
    val wrData = Bits(INPUT, C.WORDLENGTH)
    val boot   = Bool(INPUT)
    val out    = Bits(OUTPUT, C.WORDLENGTH)
    val pc    = UInt(OUTPUT, C.WORDLENGTH)
    val inst    = UInt(OUTPUT, C.WORDLENGTH)
  }

  // instantiate datapath and do instruction decoding 
  val dp = Module(new Dpath)
  dp.io.inst <> io.inst
  val zero = dp.io.zero
  def op(o: UInt) = (io.inst(31,26) === o)
  def fop(o: UInt, f: UInt) = (io.inst(31,26) === o && io.inst(5, 0) === f)

  // control logic
  val alu_src = (op(C.OP_RTYPE) || op(C.OP_BEQ) || op(C.OP_BNE))
  val alu_op = MuxCase(UInt(0), Array(
                (fop(C.OP_RTYPE, C.FUNC_ADDU) || 
                  op(C.OP_ADDI) || op(C.OP_ADDIU) || op(C.OP_SW)
                  || op(C.OP_LW)) -> (C.ALU_ADD),
                (fop(C.OP_RTYPE, C.FUNC_SUBU) || 
                  op(C.OP_BNE) || op(C.OP_BEQ)) -> (C.ALU_SUB),
                fop(C.OP_RTYPE, C.FUNC_SLTU) -> (C.ALU_SLT),
                (fop(C.OP_RTYPE, C.FUNC_OR) || op(C.OP_ORI)) -> (C.ALU_OR),
                op(C.OP_LUI) -> (C.ALU_LUIB)))
  val reg_write = (op(C.OP_RTYPE) || op(C.OP_ADDI) || op(C.OP_ADDIU)
                || op(C.OP_JAL) || op(C.OP_LW) || op(C.OP_ORI) || op(C.OP_LUI))
  val reg_dst = op(C.OP_RTYPE)
  val mem_toreg = op(C.OP_LW)
  val mem_write = op(C.OP_SW)
  val jump = (op(C.OP_J) || op(C.OP_JAL) || fop(C.OP_JR, C.FUNC_JR))
  val jrmux = fop(C.OP_JR, C.FUNC_JR)
  val branch = (op(C.OP_BEQ) || op(C.OP_BNE))
  val link_write = op(C.OP_JAL)

  // wire up the datapath 
  dp.io.isWr <> io.isWr
  dp.io.wrAddr <> io.wrAddr
  dp.io.wrData <> io.wrData
  dp.io.boot <> io.boot
  dp.io.out <> io.out
  dp.io.pc <> io.pc
  dp.io.alu_src := alu_src
  dp.io.alu_op := alu_op
  dp.io.reg_write := reg_write
  dp.io.reg_dst := reg_dst
  dp.io.mem_toreg := mem_toreg
  dp.io.mem_write := mem_write
  dp.io.jump := jump
  dp.io.jrmux := jrmux
  dp.io.branch := branch
  dp.io.link_write := link_write

}

class SingleTests(c: Single) extends Tester(c) {  
  def wr(addr: UInt, data: UInt)  = {
    poke(c.io.isWr,   1)
    poke(c.io.wrAddr, addr.litValue())
    poke(c.io.wrData, data.litValue())
    step(1)
  }
  def boot()  = {
    poke(c.io.isWr, 0)
    poke(c.io.boot, 1)
    peek(c.io.out)
    step(1)
    poke(c.io.isWr, 0)
    poke(c.io.boot, 0)
    peek(c.io.out)
    step(1)
  }
  def tick()  = {
    step(1)
    peek(c.io.out)
    peek(c.io.inst)
  }
  def IR (op: Int, rsi: Int, rti: Int, rdi: Int, shamt: Int, func: Int) = 
    Cat(UInt(op, 6), UInt(rsi, 5), UInt(rti, 5), 
        UInt(rdi, 5), UInt(shamt, 5), UInt(func, 5))
  def II (op: Int, rsi: Int, rti: Int, imm: Int) = 
    Cat(UInt(op, 6), UInt(rsi, 5), UInt(rti, 5), UInt(imm, 16))
  def StrUInt(s: String) = UInt(BigInt(s, 16))

//CHANGE HEX INST!!!//
  val app  = Array(
                    StrUInt("8fa40000"), //  lw $4, 0($29)
                    StrUInt("27a50004"), //  addiu $5, $29, 4
                    StrUInt("24a60004"), //  addiu $6, $5, 4
                    StrUInt("00041080"), //  sll $2, $4, 2
                    StrUInt("00c23021"), //  addu $6, $6, $2
                    StrUInt("0c100009"), //  jal 0x00400024 [main]	//line 0x14
                    StrUInt("00000000"), //  nop
                    StrUInt("3402000a"), //  ori $2, $0, 10
                    StrUInt("0000000c"), //  syscall
                    StrUInt("3c01ffff"), //  lui $1, -1			//0x00400024
                    StrUInt("343d2000"), //  ori $29, $1, 8192
                    StrUInt("27bdffe0"), //  addiu $29, $29, -32
                    StrUInt("afbf0014"), //  sw $31, 20($29)
                    StrUInt("afbe0010"), //  sw $30, 16($29)
                    StrUInt("27be001c"), //  addiu $30, $29, 28
                    StrUInt("34100000"), //  ori $16, $0, 0
                    StrUInt("34110014"), //  ori $17, $0, 20
                    StrUInt("00102021"), //  addu $4, $0, $16
                    StrUInt("0c10001c"), //  jal 0x00400070 [fib]	//line 0x48
                    StrUInt("0002c821"), //  addu $25, $0, $2
                    StrUInt("22100001"), //  addi $16, $16, 1
                    StrUInt("0211402b"), //  slt $8, $16, $17
                    StrUInt("1500fffa"), //  bne $8, $0, -24 [l1-0x00400058]
                    StrUInt("00051021"), //  addu $2, $0, $5
                    StrUInt("8fbf0014"), //  lw $31, 20($29)
                    StrUInt("8fbe0010"), //  lw $30, 16($29)
                    StrUInt("27bd0020"), //  addiu $29, $29, 32
                    StrUInt("03e00008"), //  jr $31
                    StrUInt("27bdffe0"), //  addiu $29, $29, -32	//0x00400070
                    StrUInt("afbf0014"), //  sw $31, 20($29)
                    StrUInt("afbe0010"), //  sw $30, 16($29)
                    StrUInt("27be001c"), //  addiu $30, $29, 28
                    StrUInt("afc40000"), //  sw $4, 0($30)
                    StrUInt("8fc20000"), //  lw $2, 0($30)
                    StrUInt("14400003"), //  bne $2, $0, 12 [$L2-0x00400088]
                    StrUInt("34020000"), //  ori $2, $0, 0
                    StrUInt("08100034"), //  j 0x004000d0 [$L1]
                    StrUInt("34080001"), //  ori $8, $0, 1
                    StrUInt("0102482b"), //  slt $9, $8, $2
                    StrUInt("15200002"), //  blwne $9, $0, 8 [$L3-0x0040009c]
                    StrUInt("08100034"), //  j 0x004000d0 [$L1]
                    StrUInt("8fc30000"), //  lw $3, 0($30)
                    StrUInt("2462ffff"), //  addiu $2, $3, -1
                    StrUInt("00022021"), //  addu $4, $0, $2
                    StrUInt("0c10001c"), //  jal 0x00400070 [fib]	//line 0xB0
                    StrUInt("afc20004"), //  sw $2, 4($30)
                    StrUInt("8fc30000"), //  lw $3, 0($30)
                    StrUInt("2462fffe"), //  addiu $2, $3, -2
                    StrUInt("00022021"), //  addu $4, $0, $2
                    StrUInt("0c10001c"), //  jal 0x00400070 [fib]
                    StrUInt("8fc30004"), //  lw $3, 4($30)
                    StrUInt("00431021"), //  addu $2, $2, $3
                    StrUInt("8fbf0014"), //  lw $31, 20($29)		//0x004000d0
                    StrUInt("8fbe0010"), //  lw $30, 16($29)
                    StrUInt("27bd0020"), //  addiu $29, $29, 32
                    StrUInt("03e00008")) //  jr $31
  wr(UInt(0), Bits(0)) // skip reset
  for (addr <- 0 until app.length) 
    wr(UInt(addr), app(addr))
  boot()
  var k = 0
  do {
    tick(); k += 1
  } while (k < 1000000)
  expect(c.io.out, 0xa18)
}

object TutorialExamples {
  def main(args: Array[String]): Unit = {
    val tutArgs = args.slice(1, args.length)
    val res =
    args(0) match {
      case "Single" =>
        chiselMainTest(tutArgs, () => Module(new Single())){
          c => new SingleTests(c)}
    }
  }
}

