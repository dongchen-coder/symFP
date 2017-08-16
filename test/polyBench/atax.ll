; ModuleID = 'polyBench/atax.bc'
source_filename = "/Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/polyBench/atax.c"
target datalayout = "e-m:o-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-apple-macosx10.12.0"

; Function Attrs: noinline nounwind ssp uwtable
define void @atax_cpu(i32 %nx, i32 %ny, double* %A, double* %x, double* %y, double* %tmp) #0 {
entry:
  %nx.addr = alloca i32, align 4
  %ny.addr = alloca i32, align 4
  %A.addr = alloca double*, align 8
  %x.addr = alloca double*, align 8
  %y.addr = alloca double*, align 8
  %tmp.addr = alloca double*, align 8
  %i = alloca i32, align 4
  %j = alloca i32, align 4
  store i32 %nx, i32* %nx.addr, align 4
  store i32 %ny, i32* %ny.addr, align 4
  store double* %A, double** %A.addr, align 8
  store double* %x, double** %x.addr, align 8
  store double* %y, double** %y.addr, align 8
  store double* %tmp, double** %tmp.addr, align 8
  store i32 0, i32* %i, align 4
  br label %for.cond

for.cond:                                         ; preds = %for.inc, %entry
  %0 = load i32, i32* %i, align 4
  %cmp = icmp slt i32 %0, 1024
  br i1 %cmp, label %for.body, label %for.end

for.body:                                         ; preds = %for.cond
  %1 = load double*, double** %y.addr, align 8
  %2 = load i32, i32* %i, align 4
  %idxprom = sext i32 %2 to i64
  %arrayidx = getelementptr inbounds double, double* %1, i64 %idxprom
  store double 0.000000e+00, double* %arrayidx, align 8
  br label %for.inc

for.inc:                                          ; preds = %for.body
  %3 = load i32, i32* %i, align 4
  %inc = add nsw i32 %3, 1
  store i32 %inc, i32* %i, align 4
  br label %for.cond

for.end:                                          ; preds = %for.cond
  store i32 0, i32* %i, align 4
  br label %for.cond3

for.cond3:                                        ; preds = %for.inc42, %for.end
  %4 = load i32, i32* %i, align 4
  %cmp4 = icmp slt i32 %4, 1024
  br i1 %cmp4, label %for.body5, label %for.end44

for.body5:                                        ; preds = %for.cond3
  %5 = load double*, double** %tmp.addr, align 8
  %6 = load i32, i32* %i, align 4
  %idxprom6 = sext i32 %6 to i64
  %arrayidx7 = getelementptr inbounds double, double* %5, i64 %idxprom6
  store double 0.000000e+00, double* %arrayidx7, align 8
  store i32 0, i32* %j, align 4
  br label %for.cond8

for.cond8:                                        ; preds = %for.inc21, %for.body5
  %7 = load i32, i32* %j, align 4
  %cmp9 = icmp slt i32 %7, 1024
  br i1 %cmp9, label %for.body10, label %for.end23

for.body10:                                       ; preds = %for.cond8
  %8 = load double*, double** %tmp.addr, align 8
  %9 = load i32, i32* %i, align 4
  %idxprom11 = sext i32 %9 to i64
  %arrayidx12 = getelementptr inbounds double, double* %8, i64 %idxprom11
  %10 = load double, double* %arrayidx12, align 8
  %11 = load double*, double** %A.addr, align 8
  %12 = load i32, i32* %i, align 4
  %mul = mul nsw i32 %12, 1024
  %13 = load i32, i32* %j, align 4
  %add = add nsw i32 %mul, %13
  %idxprom13 = sext i32 %add to i64
  %arrayidx14 = getelementptr inbounds double, double* %11, i64 %idxprom13
  %14 = load double, double* %arrayidx14, align 8
  %15 = load double*, double** %x.addr, align 8
  %16 = load i32, i32* %j, align 4
  %idxprom15 = sext i32 %16 to i64
  %arrayidx16 = getelementptr inbounds double, double* %15, i64 %idxprom15
  %17 = load double, double* %arrayidx16, align 8
  %mul17 = fmul double %14, %17
  %add18 = fadd double %10, %mul17
  %18 = load double*, double** %tmp.addr, align 8
  %19 = load i32, i32* %i, align 4
  %idxprom19 = sext i32 %19 to i64
  %arrayidx20 = getelementptr inbounds double, double* %18, i64 %idxprom19
  store double %add18, double* %arrayidx20, align 8
  br label %for.inc21

for.inc21:                                        ; preds = %for.body10
  %20 = load i32, i32* %j, align 4
  %inc22 = add nsw i32 %20, 1
  store i32 %inc22, i32* %j, align 4
  br label %for.cond8

for.end23:                                        ; preds = %for.cond8
  store i32 0, i32* %j, align 4
  br label %for.cond24

for.cond24:                                       ; preds = %for.inc39, %for.end23
  %21 = load i32, i32* %j, align 4
  %cmp25 = icmp slt i32 %21, 1024
  br i1 %cmp25, label %for.body26, label %for.end41

for.body26:                                       ; preds = %for.cond24
  %22 = load double*, double** %y.addr, align 8
  %23 = load i32, i32* %j, align 4
  %idxprom27 = sext i32 %23 to i64
  %arrayidx28 = getelementptr inbounds double, double* %22, i64 %idxprom27
  %24 = load double, double* %arrayidx28, align 8
  %25 = load double*, double** %A.addr, align 8
  %26 = load i32, i32* %i, align 4
  %mul29 = mul nsw i32 %26, 1024
  %27 = load i32, i32* %j, align 4
  %add30 = add nsw i32 %mul29, %27
  %idxprom31 = sext i32 %add30 to i64
  %arrayidx32 = getelementptr inbounds double, double* %25, i64 %idxprom31
  %28 = load double, double* %arrayidx32, align 8
  %29 = load double*, double** %tmp.addr, align 8
  %30 = load i32, i32* %i, align 4
  %idxprom33 = sext i32 %30 to i64
  %arrayidx34 = getelementptr inbounds double, double* %29, i64 %idxprom33
  %31 = load double, double* %arrayidx34, align 8
  %mul35 = fmul double %28, %31
  %add36 = fadd double %24, %mul35
  %32 = load double*, double** %y.addr, align 8
  %33 = load i32, i32* %j, align 4
  %idxprom37 = sext i32 %33 to i64
  %arrayidx38 = getelementptr inbounds double, double* %32, i64 %idxprom37
  store double %add36, double* %arrayidx38, align 8
  br label %for.inc39

for.inc39:                                        ; preds = %for.body26
  %34 = load i32, i32* %j, align 4
  %inc40 = add nsw i32 %34, 1
  store i32 %inc40, i32* %j, align 4
  br label %for.cond24

for.end41:                                        ; preds = %for.cond24
  br label %for.inc42

for.inc42:                                        ; preds = %for.end41
  %35 = load i32, i32* %i, align 4
  %inc43 = add nsw i32 %35, 1
  store i32 %inc43, i32* %i, align 4
  br label %for.cond3

for.end44:                                        ; preds = %for.cond3
  ret void
}

attributes #0 = { noinline nounwind ssp uwtable "correctly-rounded-divide-sqrt-fp-math"="false" "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-jump-tables"="false" "no-nans-fp-math"="false" "no-signed-zeros-fp-math"="false" "no-trapping-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="false" "use-soft-float"="false" }

!llvm.module.flags = !{!0}
!llvm.ident = !{!1}

!0 = !{i32 1, !"PIC Level", i32 2}
!1 = !{!"clang version 4.0.0 (tags/RELEASE_400/final)"}
