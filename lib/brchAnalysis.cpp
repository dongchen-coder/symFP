//
//  brchAnalysis.cpp
//  LLVM
//
//  Created by Fangzhou Liu on 2017/7/11.
//
//


#include "brchAnalysis.hpp"

using namespace std;
using namespace llvm;


namespace brchAnalysis {
    bool isLoopCondition(Loop* L, BasicBlock* bb);
    char BranchAnalysis::ID = 0;
    static RegisterPass<BranchAnalysis> X("brchAnalysis", "branch condition analysis Pass", false, false);
    
    BranchAnalysis::BranchAnalysis() : FunctionPass(ID) {}
    
    struct LoopBranchStruct {
        Loop *L;
        vector<ICmpInst *> BR;
        vector<Value *> ARR;
    };
    
    vector<LoopBranchStruct> LoopBranchInfo;
    
    
    vector<BasicBlock *> BranchAnalysis::FindIfBody(Loop *L) {
        vector<BasicBlock *> body;
        bool begin = false;
        for (BasicBlock *BB: L->getBlocks()) {
            //errs() << "BB: " << BB->getName() << "\n";
            //errs() << "BB: " << BB->getName() << "\n";
            if (std::regex_match (BB->getName().str(), std::regex("^if.then$|^if.then\\d*$)"))) {
                begin = true;
            }
            if (begin) {
                body.push_back(BB);
            }
            if (std::regex_match (BB->getName().str(), std::regex("^if.end$|^if.end\\d*$)"))) {
                begin = false;
            }
        }
        return body;
    }
    
    bool BranchAnalysis::runOnFunction(Function &F) {
        
        errs() << "\nStart analysis Branch Conditions \n";
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* lrn = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;

//        getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().DumpLoopTree(lrn, "");

        // tranverse the LoopRefNode
        // The root of the LoopRefTree is a special node with all component except the 'next' is empty
        if (lrn != NULL && lrn->next != NULL) {
            for (TreeNodeBase* node: *(lrn->next)) {
                tranverseLoopRefTree(node);
            }
        }
        

        // Can not iterpret the info
//        dumpBranchInfoStruct();
        
        return false;
    }

    void BranchAnalysis::tranverseLoopRefTree(TreeNodeBase* node) {
        if (node != NULL) {
            if (static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(node) == nullptr) {
                return;
            }
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * lnd = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(node);
            if (lnd->L != NULL) {
                errs() << lnd->L->getName();
                FindBranch(lnd);
                vector<TreeNodeBase *>* subnodes = lnd->next;
                for (vector<TreeNodeBase *>::iterator nit=subnodes->begin(), neit=subnodes->end(); nit!=neit; ++nit) {
                    tranverseLoopRefTree(*nit);
                }
            }
        }
    }

    void BranchAnalysis::FindBranch(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node) {
        Loop* loop = node->L;
        loopAnalysis::LoopIndvBoundAnalysis::LoopInfoStruct* lstr = node->LIS;
//        errs() << "\nFind branch conditions in Loop " << loop->getName() << "\n";
        for (BasicBlock* loopbb: loop->getBlocks()) {
            for (BasicBlock::iterator it=loopbb->begin(), eit=loopbb->end(); it!=eit; ++it) {
                if (isa<BranchInst>(*it)) {
                    BranchInst* brinst = dyn_cast<BranchInst>(it);
                    if (brinst->isConditional()) {
                        if  (isLoopCondition(loop, loopbb)) {
//                            errs() << "Is Loop Condition \n";
                        }
                        else {
                            Value *cond = brinst->getCondition();
//                            errs() << "Condition is " << *cond << " => " << dumpBranchCondition(cond) << "\n";
                            if (isa<ICmpInst>(cond)) {
                                ICmpInst* cmpinst = dyn_cast<ICmpInst>(cond);
                                // Create the condition info struct
                                ConditionInfoStruct temp = {
                                    node,
                                    loop,
                                    loopbb,
                                    cmpinst
//                                    cmpinst->getOperand(0),
//                                    cmpinst->getOperand(1),
//                                    cmpinst->getPredicate()
                                };
                                dumpCondInfoStruch(temp);
                            }
                        }
                    }
                }
            }
        }
//        vector<Value *> arrlist;
//        vector<ICmpInst *> condlist;
//        vector<BasicBlock*> body = FindIfBody(L);
//        for (BasicBlock *b: body) {
//            for(BasicBlock::iterator II = b->begin(); II != b->end(); ++II) {
//                // has array visit
//                if (isa<GetElementPtrInst>(*II)) {
//                    arrlist.push_back(II->getOperand(0));
//                }
//            }
//        }
//        // if there's array visit in the if-body
//        if (arrlist.size() > 0) {
//            BasicBlock *begin = body.front();
//            // find the if condition
//            // iterate each predecessor of the first BB in if-body
//            for(auto pit = pred_begin(begin), pet = pred_end(begin); pit != pet; ++pit) {
//                BasicBlock *predecessor = *pit;
//                for (BasicBlock::iterator II = predecessor->begin(); II != predecessor->end(); ++II) {
//                    if (isa<ICmpInst>(*II) && (CheckBranchValid(L, II->getOperand(0)) || CheckBranchValid(L, II->getOperand(1)))) {
//                        // condition
//                        condlist.push_back(cast<ICmpInst>(&*II));
//                    }
//                }
//            }
//            // store the branch value
//            if (condlist.size() > 0) {
//                LoopBranchStruct temp = {
//                    L,
//                    condlist,
//                    arrlist
//                };
//                LoopBranchInfo.push_back(temp);
//            }
//        }
        
    }
    
    /*
     * Check Whether the branch condition is the loop condition
     * This means that this branch condition cannot be the loop condition, or its subloop condition
     * @para *L        Loop Object
     * @para *bb       BasicBlock that contains a condition instruction
     */
    bool isLoopCondition(Loop* L, BasicBlock *bb) {
//        errs() << "BasicBlock " << bb->getName() << " in Loop " << dyn_cast<BasicBlock>(*(L->block_begin()))->getName() << "\n";
        bool result = false;
        if (bb->getName().str().compare(dyn_cast<BasicBlock>(*(L->block_begin()))->getName().str()) != 0) {
            for (Loop *subloop: L->getSubLoopsVector()) {
                result = result || isLoopCondition(subloop, bb);
            }
        } else {
            result = true;
        }
        return result;
    }
  
#warning Doesn't use
    /*
     * Check Whether the Branch Condition contains the Loop Indv
     * @para *L     Loop Object to be analyzed
     * @para val    Variables contains in the
     */
    bool BranchAnalysis::CheckBranchValid(Loop* L, Value* val) {
        if (isa<ConstantInt>(val)) { // constant value
            return false;
        }
        else if (isa<LoadInst>(val)) {
            LoadInst* load = dyn_cast<LoadInst>(val);
            return CheckIndvVar(L, load->llvm::User::getOperand(0)->getName().str());
        }
        else if (isa<Instruction>(val)) {
            Instruction *instr = dyn_cast<Instruction>(val);
            if (isa<GetElementPtrInst>(instr)) { // array in the condition
                return false;
            }
            switch (instr->getNumOperands()) {
                case 0:
                    return false;
                    break;
                case 1:
                    return CheckBranchValid(L, instr->getOperand(0));
                    break;
                case 2:
                    return (CheckBranchValid(L, instr->getOperand(0)) || CheckBranchValid(L, instr->getOperand(1)));
                    break;
                default:
                    return false;
            }
        }
        else {
            return false;
        }
    }

#warning Doesn't use
    /*
     * Check Whether the Condition Var belongs to the Loop Indv of the current Loop/its subLoop 
     * @para *L     Loop Object to be analyzed
     * @para var    Condition Variable
     */
    bool BranchAnalysis::CheckIndvVar(Loop* L, string var) {
        loopAnalysis::LoopIndvBoundAnalysis::LoopInfoStruct* ls = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ExtractLoopInfo(L);

        // check whether the condition var belongs to loop IDV of the current loop
        for (vector<Value* >::iterator it=ls->IDV->begin(), eit=ls->IDV->end(); it != eit; ++it) {
            Value* idv = dyn_cast<Value>(*it);
            if (idv->getName().str().compare(var) == 0) {
                return true;
            }
        }
    return false;
    }

    
    void BranchAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        return;
    }

#warning Doesn't use
    void BranchAnalysis::dumpBranchInfoStruct() {
        if (LoopBranchInfo.size() == 0) {
            errs() << "\n ====== No Satisfied Loop Branch Conditions ======\n";
            return;
        }
        for (LoopBranchStruct lbs: LoopBranchInfo) {
            errs() << "\n ====== For Loop " << lbs.L->getName() << " ======\n";
            errs() << "Condition Branch: \t";
            for (ICmpInst *br: lbs.BR) {
                switch (br->getPredicate()) {
                    case ICmpInst::ICMP_EQ:
                        errs() << "(" << dumpValue(br->getOperand(0)) << " == " << dumpValue(br->getOperand(1)) << ")";
                        break;
                    case ICmpInst::ICMP_NE:
                        errs() << "(" << dumpValue(br->getOperand(0)) << " != " << dumpValue(br->getOperand(1)) << ")";
                        break;
                    case ICmpInst::ICMP_SGE:
                        errs() << "(" << dumpValue(br->getOperand(0)) << " >= " << dumpValue(br->getOperand(1)) << ")";
                        break;
                    case ICmpInst::ICMP_SGT:
                        errs() << "(" << dumpValue(br->getOperand(0)) << " > " << dumpValue(br->getOperand(1)) << ")";
                        break;
                    case ICmpInst::ICMP_SLE:
                        errs() << "(" << dumpValue(br->getOperand(0)) << " <= " << dumpValue(br->getOperand(1)) << ")";
                        break;
                    case ICmpInst::ICMP_SLT:
                        errs() << "(" << dumpValue(br->getOperand(0)) << " < " << dumpValue(br->getOperand(1)) << ")";
                        break;
                    default:
                        break;
                }
                errs() << " ";
            }
            errs() << "\nBranch Arrays: \t";
            for (Value *arr: lbs.ARR) {
                errs() << dumpValue(arr) << " ";
            }
            errs() << "\n ====== End ======\n\n";
        }

        
    }
    
    string BranchAnalysis::dumpBranchCondition(Value* cond) {
        if (isa<ICmpInst>(cond)) {
            ICmpInst *br = dyn_cast<ICmpInst>(cond);
            
            switch (br->getPredicate()) {
                case ICmpInst::ICMP_EQ:
                    return "(" + dumpValue(br->getOperand(0)) + " == " + dumpValue(br->getOperand(1)) + ")";
                    break;
                case ICmpInst::ICMP_NE:
                    return "(" + dumpValue(br->getOperand(0)) + " != " + dumpValue(br->getOperand(1)) + ")";
                    break;
                case ICmpInst::ICMP_SGE:
                    return "(" + dumpValue(br->getOperand(0)) + " >= " + dumpValue(br->getOperand(1)) + ")";
                    break;
                case ICmpInst::ICMP_SGT:
                    return  "(" + dumpValue(br->getOperand(0)) + " > " + dumpValue(br->getOperand(1)) + ")";
                    break;
                case ICmpInst::ICMP_SLE:
                    return "(" + dumpValue(br->getOperand(0)) + " <= " + dumpValue(br->getOperand(1)) + ")";
                    break;
                case ICmpInst::ICMP_SLT:
                    return "(" + dumpValue(br->getOperand(0)) + " < " + dumpValue(br->getOperand(1)) + ")";
                    break;
                default:
                    break;
            }
        }
        return dumpValue(cond);
    }

    
    string BranchAnalysis::dumpValue(Value *v) {
        
        if (isa<Instruction>(v)) {
            
            Instruction *inst = cast<Instruction>(v);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    return "(" + dumpValue(inst->getOperand(0)) + " + " + dumpValue(inst->getOperand(1)) + ")";
                    break;
                case Instruction::Sub:
                    return "(" + dumpValue(inst->getOperand(0)) + " - " + dumpValue(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Mul:
                    return "(" + dumpValue(inst->getOperand(0)) + " * " + dumpValue(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return "(" + dumpValue(inst->getOperand(0)) + " / " + dumpValue(inst->getOperand(1)) + ")";
                    break;
                case Instruction::Load:
                    return cast<LoadInst>(inst)->getPointerOperand()->getName().str();
                    break;
                case Instruction::GetElementPtr:
                    return dumpArray(inst->getOperand(0));
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(v)) {
            return to_string(dyn_cast<ConstantInt>(v)->getValue().getSExtValue());
        }
        return "";

        
    }
    
#warning Doesn't use
    string BranchAnalysis::dumpArray(Value *arr) {
        if (isa<Instruction>(arr)) {
            Instruction *arrayInst = cast<Instruction>(arr);
            switch (arrayInst->getOpcode()) {
                case Instruction::SExt:
                    // this is the array index
                    break;
                case Instruction::Load:
                    return cast<LoadInst>(arrayInst)->getPointerOperand()->getName().str();
                    break;
                default:
                    break;
            }
        }
        return "";
    }
    
    void BranchAnalysis::dumpCondInfoStruch(ConditionInfoStruct cis) {
        errs() << "===========================================\n";
        errs() << "In Loop: \t" << cis.L->getName() << "\n";
        errs() << "In BasicBlock: \t" << cis.BB->getName() << "\n";
        errs() << "Condition: \t" << dumpBranchCondition(cis.COND) << "\n";
        errs() << "===========================================\n";
    }
}

